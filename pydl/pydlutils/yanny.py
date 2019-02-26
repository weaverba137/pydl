(str,)# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the yanny directory in idlutils.

This is a Python library for reading & writing yanny files.

:class:`yanny` is an object-oriented interface to SDSS Parameter files
(a.k.a. "FTCL" or "yanny" files) following these specifications_.
Parameter files typically have and can be recognized by the extension
``.par``.  These files may also be recognized if the first line of the
file is::

    #%yanny

This is not part of the standard specification, but it suffices to
identify, *e.g.*, files that have not yet been written to disk, but only
exist as file objects.

Because Parameter files can contain multiple tables, as well as
metadata, there is no simple, one-to-one correspondence between these
files and, say, an astropy :class:`~astropy.table.Table` object.  Thus,
:class:`yanny` objects are implemented as a subclass of
:class:`~collections.OrderedDict` (to remember the order of keyword-value
pairs), and the actual data are values accessed by keyword.  Still, it is
certainly possible to *write* a table-like object to a yanny file.

Given the caveats above, we have introduced *experimental* support for
reading and writing yanny files directly to/from :class:`~astropy.table.Table`
objects.  Because of the experimental nature of this support, it must be
activated "by hand" by including this snippet in your code::

    from astropy.table import Table
    from astropy.io.registry import (register_identifier, register_reader,
                                     register_writer)
    from pydl.pydlutils.yanny import (is_yanny, read_table_yanny,
                                      write_table_yanny)

    register_identifier('yanny', Table, is_yanny)
    register_reader('yanny', Table, read_table_yanny)
    register_writer('yanny', Table, write_table_yanny)

Currently multidimensional arrays are only supported for type ``char``, and a
close reading of the specifications indicates that multidimensional arrays
were only ever intended to be supported for type ``char``.  So no
multidimensional arrays, sorry.

.. _specifications: https://www.sdss.org/dr14/software/par/
"""
import re
import os
import datetime
import warnings
from collections import OrderedDict
import numpy as np
from astropy.table import Table
# from astropy.io.registry import register_identifier, register_writer
from . import PydlutilsException, PydlutilsUserWarning


class yanny(OrderedDict):
    """An object interface to a yanny file.

    Create a yanny object using a yanny file, `filename`.  If the file exists,
    it is read, & the dict structure of the object will be basically the
    same as that returned by ``read_yanny()`` in the efftickle package.

    If the file does not exist, or if no filename is given, a blank
    structure is returned.  Other methods allow for subsequent writing
    to the file.

    Parameters
    ----------
    filename : :class:`str` or file-like, optional
        The name of a yanny file or a file-like object representing a
        yanny file.
    raw : :class:`bool`, optional
        If ``True``, data in a yanny file will *not* be converted into
        astropy :class:`~astropy.table.Table` objects, but will instead be
        retained as raw Python lists.

    Attributes
    ----------
    raw : :class:`bool`
        If ``True``, data in a yanny file will *not* be converted into
        astropy :class:`~astropy.table.Table` objects, but will instead be
        retained as raw Python lists.
    filename : :class:`str`
        The name of a yanny parameter file.  If a file-like object was used
        to initialize the object, this will have the value 'in_memory.par'.
    _symbols : :class:`dict`
        A dictionary containing the metadata describing the tables.
    _contents : :class:`str`
        The complete contents of a yanny parameter file.
    _struct_type_caches : :class:`dict`
        A dictionary of dictionaries, one dictionary for every structure
        definition in a yanny parameter file.  Contains the types of
        each column
    _struct_isarray_caches : :class:`dict`
        A dictionary of dictionaries, one dictionary for every structure
        definition in a yanny parameter file.  Contains a boolean value
        for every column.
    _enum_cache : :class:`dict`
        Initially ``None``, this attribute is initialized the first time
        the :meth:`isenum` method is called.  The keyword is the name of the
        enum type, the value is a list of the possible values of that type.

    """

    @staticmethod
    def get_token(string):
        """Removes the first 'word' from string.

        If the 'word' is enclosed in double quotes, it returns the
        contents of the double quotes. If the 'word' is enclosed in
        braces, it returns the contents of the braces, but does not
        attempt to split the array.  If the 'word' is the last word of the
        string, remainder is set equal to the empty string.  This is
        basically a wrapper on some convenient regular expressions.

        Parameters
        ----------
        string : :class:`str`
            A string containing words.

        Returns
        -------
        :func:`tuple`
            A tuple containing the first word and the remainder of the string.

        Examples
        --------
        >>> from pydl.pydlutils.yanny import yanny
        >>> yanny.get_token("The quick brown fox")
        ('The', 'quick brown fox')
        """
        if string[0] == '"':
            (word, remainder) = re.search(r'^"([^"]*)"\s*(.*)',
                                          string).groups()
        elif string[0] == '{':
            (word, remainder) = re.search(r'^\{\s*([^}]*)\s*\}\s*(.*)',
                                          string).groups()
        else:
            try:
                (word, remainder) = re.split(r'\s+', string, 1)
            except ValueError:
                (word, remainder) = (string, '')
        # if remainder is None:
        #     remainder = ''
        return (word, remainder)

    @staticmethod
    def protect(x):
        """Used to appropriately quote string that might contain whitespace.

        This method is mostly for internal use by the yanny object.

        Parameters
        ----------
        x : :class:`str`
            The data to protect.

        Returns
        -------
        :class:`str`
            The data with white space protected by quotes.

        Examples
        --------
        >>> from pydl.pydlutils.yanny import yanny
        >>> yanny.protect('This string contains whitespace.')
        '"This string contains whitespace."'
        >>> yanny.protect('This string contains a #hashtag.')
        '"This string contains a #hashtag."'
        """
        if isinstance(x, np.bytes_):
            s = x.decode()
        else:
            s = str(x)
        if len(s) == 0 or s.find('#') >= 0 or re.search(r'\s+', s) is not None:
            return '"' + s + '"'
        else:
            return s

    @staticmethod
    def trailing_comment(line):
        """Identify a trailing comment and strip it.

        This routine works on the theory that a properly quoted comment mark
        will be surrounted by an odd number of double quotes, & we can
        skip to searching for the last one in the line.

        Parameters
        ----------
        line : :class:`str`
            A line from a yanny file potentially containing trailing comments.

        Returns
        -------
        :class:`str`
            The line with any trailing comment and any residual white space
            trimmed off.

        Notes
        -----
        This may fail in certain pathological cases, for example if a
        real trailing comment contains a single double-quote::

            # a 'pathological" trailing comment

        or if someone is over-enthusiastically commenting::

            # # # # # I like # characters.

        Examples
        --------
        >>> from pydl.pydlutils.yanny import yanny
        >>> yanny.trailing_comment('mystruct 1234 "#hashtag" # a comment.')
        'mystruct 1234 "#hashtag"'
        >>> yanny.trailing_comment('mystruct 1234 "#hashtag" # a "comment".')
        'mystruct 1234 "#hashtag"'
        """
        lastmark = line.rfind('#')
        if lastmark >= 0:
            #
            # Count the number of double quotes in the remainder of the line
            #
            if (len([c for c in line[lastmark:] if c == '"']) % 2) == 0:
                #
                # Even number of quotes
                #
                return line[0:lastmark].rstrip()
        return line

    @staticmethod
    def dtype_to_struct(dt, structname='mystruct', enums=None):
        """Convert a NumPy dtype object describing a record array to
        a typedef struct statement.

        The second argument is the name of the structure.
        If any of the columns are enum types, enums must
        be a dictionary with the keys the column names, and the values
        are a tuple containing the name of the enum type as the first item
        and a tuple or list of possible values as the second item.

        Parameters
        ----------
        dt : :class:`numpy.dtype`
            The dtype of a NumPy record array.
        structname : :class:`str`, optional
            The name to give the structure in the yanny file.  Defaults to
            'MYSTRUCT'.
        enums : :class:`dict`, optional
            A dictionary containing enum information.  See details above.

        Returns
        -------
        :class:`dict`
            A dictionary suitable for setting the 'symbols' dictionary of a
            new yanny object.

        Examples
        --------
        """
        dtmap = {'i2': 'short', 'i4': 'int', 'i8': 'long', 'f4': 'float',
                 'f8': 'double'}
        returnenums = list()
        if enums is None:
            enums = dict()
        else:
            for e in enums:
                lines = list()
                lines.append('typedef enum {')
                for n in enums[e][1]:
                    lines.append("    {0},".format(n))
                lines[-1] = lines[-1].strip(',')
                lines.append('}} {0};'.format(enums[e][0].upper()))
                returnenums.append("\n".join(lines))
                # lines.append('')
        lines = list()
        lines.append('typedef struct {')
        for c in dt.names:
            if dt[c].kind == 'V':
                t = dt[c].subdtype[0].str[1:]
                l = dt[c].subdtype[1][0]
                s = dt[c].subdtype[0].itemsize
            else:
                t = dt[c].str[1:]
                l = 0
                s = dt[c].itemsize
            line = '    '
            if t[0] in 'SU':
                if c in enums:
                    line += enums[c][0].upper()
                else:
                    line += 'char'
            else:
                line += dtmap[t]
            line += ' {0}'.format(c)
            if l > 0:
                line += "[{0:d}]".format(l)
            if t[0] in 'SU' and c not in enums:
                line += "[{0:d}]".format(s)
            line += ';'
            lines.append(line)
        lines.append('}} {0};'.format(structname.upper()))
        return {structname.upper(): list(dt.names),
                'enum': returnenums, 'struct': ["\n".join(lines)]}

    def __init__(self, filename=None, raw=False):
        """Create a yanny object using a yanny file.
        """
        super(yanny, self).__init__()
        #
        # The symbol hash is inherited from the old read_yanny
        #
        self._symbols = dict()
        #
        # Create special attributes that contain the internal status of the
        # object. This should prevent overlap with keywords in the data files.
        #
        self.filename = ''
        self._contents = ''
        #
        # Since the re is expensive, cache the structure types keyed by
        # the field. Create a dictionary for each structure found.
        #
        self._struct_type_caches = dict()
        self._struct_isarray_caches = dict()
        self._enum_cache = None
        #
        # Optionally convert numeric data into NumPy arrays
        #
        self.raw = raw
        #
        # If the file exists, read it
        #
        if filename is not None:
            #
            # Handle file-like objects
            #
            if isinstance(filename, (str,)):
                if os.access(filename, os.R_OK):
                    self.filename = filename
                    with open(filename, 'r') as f:
                        self._contents = f.read()
            else:
                #
                # Assume file-like
                #
                self.filename = 'in_memory.par'
                contents = filename.read()
                if 'b' in filename.mode:
                    contents = contents.decode('ascii')
                self._contents = contents
            self._parse()
        return

    def __str__(self):
        """Implement the ``str()`` function for yanny objects.

        Simply prints the current contents of the yanny file.
        """
        return self._contents

    __repr__ = __str__

    def __eq__(self, other):
        """Test two yanny objects for equality.

        Two yanny objects are assumed to be equal if their contents are equal.
        """
        if isinstance(other, yanny):
            return self._contents == other._contents
        return NotImplemented

    def __ne__(self, other):
        """Test two yanny objects for inequality.

        Two yanny objects are assumed to be unequal if their contents are
        unequal.
        """
        if isinstance(other, yanny):
            return self._contents != other._contents
        return NotImplemented

    def __bool__(self):
        """Give a yanny object a definite truth value.

        A yanny object is considered ``True`` if its contents are non-zero.
        """
        return len(self._contents) > 0

    # `__nonzero__` is needed for Python 2.
    # Python 3 uses `__bool__`.
    # http://stackoverflow.com/a/2233850/498873
    __nonzero__ = __bool__

    def type(self, structure, variable):
        """Returns the type of a variable defined in a structure.

        Returns ``None`` if the structure or the variable is undefined.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column whose type you want.

        Returns
        -------
        :class:`str`
            The type of the variable.
        """
        if structure not in self:
            return None
        if variable not in self.columns(structure):
            return None
        #
        # Added code to cache values to speed up parsing large files.
        # 2009.05.11 / Demitri Muna, NYU
        # Find (or create) the cache for this structure.
        #
        try:
            cache = self._struct_type_caches[structure]
        except KeyError:
            self._struct_type_caches[structure] = dict()
            # cache for one struct type
            cache = self._struct_type_caches[structure]
        #
        # Lookup (or create) the value for this variable
        #
        try:
            var_type = cache[variable]
        except KeyError:
            defl = [x for x in self._symbols['struct']
                    if x.find(structure.lower()) > 0]
            defu = [x for x in self._symbols['struct']
                    if x.find(structure.upper()) > 0]
            if len(defl) != 1 and len(defu) != 1:
                return None
            elif len(defl) == 1:
                definition = defl
            else:
                definition = defu
            typere = re.compile(
                r'(\S+)\s+{0}([\[<].*[\]>]|);'.format(variable))
            (typ, array) = typere.search(definition[0]).groups()
            var_type = typ + array.replace('<', '[').replace('>', ']')
            cache[variable] = var_type
        return var_type

    def basetype(self, structure, variable):
        """Returns the bare type of a variable, stripping off any array
        information.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column whose type you want.

        Returns
        -------
        :class:`str`
            The type of the variable, stripped of array information.
        """
        typ = self.type(structure, variable)
        try:
            return typ[0:typ.index('[')]
        except ValueError:
            return typ

    def isarray(self, structure, variable):
        """Returns ``True`` if the variable is an array type.

        For character types, this means a two-dimensional array,
        *e.g.*: ``char[5][20]``.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column to check for array type.

        Returns
        -------
        :class:`bool`
            ``True`` if the variable is an array.
        """
        try:
            cache = self._struct_isarray_caches[structure]
        except KeyError:
            self._struct_isarray_caches[structure] = dict()
            cache = self._struct_isarray_caches[structure]
        try:
            result = cache[variable]
        except KeyError:
            typ = self.type(structure, variable)
            character_array = re.compile(r'char[\[<]\d*[\]>][\[<]\d*[\]>]')
            if ((character_array.search(typ) is not None) or
                    (typ.find('char') < 0 and (typ.find('[') >= 0 or
                                               typ.find('<') >= 0))):
                cache[variable] = True
            else:
                cache[variable] = False
            result = cache[variable]
        return result

    def isenum(self, structure, variable):
        """Returns true if a variable is an enum type.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column to check for enum type.

        Returns
        -------
        :class:`bool`
            ``True`` if the variable is enum type.
        """
        if self._enum_cache is None:
            self._enum_cache = dict()
            if 'enum' in self._symbols:
                for e in self._symbols['enum']:
                    m = re.search(r'typedef\s+enum\s*\{([^}]+)\}\s*(\w+)\s*;',
                                  e).groups()
                    self._enum_cache[m[1]] = re.split(r',\s*', m[0].strip())
            else:
                return False
        return self.basetype(structure, variable) in self._enum_cache

    def array_length(self, structure, variable):
        """Returns the length of an array type or 1 if the variable is not an
        array.

        For character types, this is the length of a two-dimensional
        array, *e.g.*, ``char[5][20]`` has length 5.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column to check for array length.

        Returns
        -------
        :class:`int`
            The length of the array variable
        """
        if self.isarray(structure, variable):
            typ = self.type(structure, variable)
            return int(typ[typ.index('[')+1:typ.index(']')])
        else:
            return 1

    def char_length(self, structure, variable):
        """Returns the length of a character field.

        *e.g.* ``char[5][20]`` is an array of 5 strings of length 20.
        Returns ``None`` if the variable is not a character type. If the
        length is not specified, *i.e.* ``char[]``, it returns the length of
        the largest string.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column to check for char length.

        Returns
        -------
        :class:`int` or None
            The length of the char variable.
        """
        typ = self.type(structure, variable)
        if typ.find('char') < 0:
            return None
        try:
            return int(typ[typ.rfind('[')+1:typ.rfind(']')])
        except ValueError:
            if self.isarray(structure, variable):
                return max([max([len(x) for x in r])
                            for r in self[structure][variable]])
            else:
                return max([len(x) for x in self[structure][variable]])

    def dtype(self, structure):
        """Returns a NumPy dtype object suitable for describing a table as a
        record array.

        Treats enums as string, which is what the IDL reader does.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure.

        Returns
        -------
        :class:`numpy.dtype`
            A dtype object suitable for describing the yanny structure as a
            record array.
        """
        dt = list()
        dtmap = {'short': 'i2', 'int': 'i4', 'long': 'i8', 'float': 'f',
                 'double': 'd'}
        for c in self.columns(structure):
            typ = self.basetype(structure, c)
            if typ == 'char':
                d = "S{0:d}".format(self.char_length(structure, c))
            elif self.isenum(structure, c):
                d = "S{0:d}".format(max([len(x) for x in
                                         self._enum_cache[typ]]))
            else:
                d = dtmap[typ]
            if self.isarray(structure, c):
                dt.append((str(c), str(d), (self.array_length(structure, c),)))
            else:
                dt.append((str(c), str(d)))
        dt = np.dtype(dt)
        return dt

    def convert(self, structure, variable, value):
        """Converts value into the appropriate (Python) type.

        * ``short`` & ``int`` are converted to Python :class:`int`.
        * ``long`` is converted to Python :class:`long`.
        * ``float`` & ``double`` are converted to Python :class:`float`.
        * Other types are not altered.

        There may be further conversions into NumPy types, but this is the
        first stage.

        Parameters
        ----------
        structure : :class:`str`
            The name of the structure that contains `variable`.
        variable : :class:`str`
            The name of the column undergoing conversion.
        value : :class:`str`
            The value contained in a particular row of `variable`.

        Returns
        -------
        :class:`int`, :class:`long`, :class:`float` or :class:`str`
            `value` converted to a Python numerical type.
        """
        intTypes = set(['short', 'int', 'long'])
        floatTypes = set(['float', 'double'])
        typ = self.basetype(structure, variable)
        if typ in intTypes:
            if self.isarray(structure, variable):
                return [int(v) for v in value]
            else:
                return int(value)
        if typ in floatTypes:
            if self.isarray(structure, variable):
                return [float(v) for v in value]
            else:
                return float(value)
        return value

    def tables(self):
        """Returns a list of all the defined structures.

        This is just the list of keys of the object with the 'internal'
        keys removed.
        """
        foo = list()
        for k in self._symbols.keys():
            if k not in ('struct', 'enum'):
                foo.append(k)
        return foo

    def columns(self, table):
        """Returns an ordered list of column names associated with a
        particular table.

        The order is the same order as they are defined in the yanny file.

        Parameters
        ----------
        table : :class:`str`
            The table whose columns are desired.

        Returns
        -------
        :class:`list`
            The list of column names.
        """
        foo = list()
        if table in self._symbols:
            return self._symbols[table]
        return foo

    def size(self, table):
        """Returns the number of rows in a table.

        Parameters
        ----------
        table : :class:`str`
            The table whose size desired.

        Returns
        -------
        :class:`int`
            The number of rows in `table`.
        """
        foo = self.columns(table)
        return len(self[table][foo[0]])

    def pairs(self):
        """Returns a list of keys to keyword/value pairs.

        Equivalent to doing ``self.keys()``, but with all the data tables &
        other control structures stripped out.
        """
        p = list()
        foo = self.tables()
        for k in self.keys():
            if k not in foo:
                p.append(k)
        return p

    def row(self, table, index):
        """Returns a list containing a single row from a specified table in
        column order.

        If index is out of range, it returns an empty list.

        If the yanny object instance is set up for NumPy record arrays, then
        a single row can be obtained with::

            row0 = par['TABLE'][0]

        Parameters
        ----------
        table : :class:`str`
            The table whose row is desired.
        index : :class:`int`
            The number of the row to return.

        Returns
        -------
        :class:`list`
            A row from `table`.
        """
        datarow = list()
        if table in self and index >= 0 and index < self.size(table):
            for c in self.columns(table):
                datarow.append(self[table][c][index])
        return datarow

    def list_of_dicts(self, table):
        """Construct a list of dictionaries.

        Takes a table from the yanny object and constructs a list object
        containing one row per entry. Each item in the list is a dictionary
        keyed by the struct value names.

        If the yanny object instance is set up for NumPy record arrays, then
        the same functionality can be obtained with::

            foo = par['TABLE'][0]['column']

        Parameters
        ----------
        table : :class:`str`
            The table to convert

        Returns
        -------
        :class:`list`
            A list containing the rows of `table` converted to :class:`dict`.
        """
        return_list = list()
        d = dict()
        # I'm assuming these are in order...
        struct_fields = self.columns(table)
        for i in range(self.size(table)):
            one_row = self.row(table, i)  # one row as a list
            j = 0
            for key in struct_fields:
                d[key] = one_row[j]
                j = j + 1
            return_list.append(dict(d))  # append a new dict (copy of d)
        return return_list

    def new_dict_from_pairs(self):
        """Returns a new dictionary of keyword/value pairs.

        The new dictionary (*i.e.*, not a yanny object) contains the keys
        that :meth:`pairs` returns. There are two reasons this is convenient:

        * the key 'symbols' that is part of the yanny object will not be
          present
        * a simple yanny file can be read with no further processing

        Returns
        -------
        :class:`~collections.OrderedDict`
            A dictionary of the keyword-value pairs that remembers the order
            in which they were defined in the file.

        Examples
        --------

        Read a yanny file and return only the pairs::

            >>> from os.path import dirname
            >>> from pydl.pydlutils.yanny import yanny
            >>> new_dict = yanny(dirname(__file__)+'/tests/t/test.par').new_dict_from_pairs()
            >>> new_dict['mjd']
            '54579'
            >>> new_dict['alpha']
            'beta gamma delta'

        added: Demitri Muna, NYU 2009-04-28
        """
        new_dictionary = OrderedDict()
        for key in self.pairs():
            new_dictionary[key] = self[key]
        return new_dictionary

    def write(self, newfile=None, comments=None):
        """Write a yanny object to a file.

        This assumes that the filename used to create the object was not that
        of a pre-existing file.  If a file of the same name is detected,
        this method will *not* attempt to overwrite it, but will print a warning.
        This also assumes that the special 'symbols' key has been properly
        created.  This will not necessarily make the file very human-readable,
        especially if the data lines are long.  If the name of a new file is
        given, it will write to the new file (assuming it doesn't exist).
        If the writing is successful, the data in the object will be updated.

        Parameters
        ----------
        newfile : :class:`str`, optional
            The name of the file to write.
        comments : :class:`str` or :class:`list` of :class:`str`, optional
            Comments that will be placed at the head of the file.  If a
            single string is passed, it will be written out verbatim, although
            a '#' character will be added if it does not already have one.
            If a list of strings is passed, comment characters will be added
            and the strings will be joined together.
        """
        if newfile is None:
            if len(self.filename) > 0:
                newfile = self.filename
            else:
                raise ValueError("No filename specified!")
        if os.access(newfile, os.F_OK):
            raise PydlutilsException(
                  "{0} exists, aborting write!".format(newfile))
        if comments is None:
            basefile = os.path.basename(newfile)
            timestamp = datetime.datetime.utcnow().strftime(
                        '%Y-%m-%d %H:%M:%S UTC')
            comments = "#\n# {0}\n#\n# Created by pydl.pydlutils.yanny.yanny\n#\n# {1}\n#\n".format(basefile, timestamp)
        else:
            if isinstance(comments, (str,)):
                if not comments.startswith('#'):
                    comments = '# ' + comments
                if not comments.endswith('\n'):
                    comments += '\n'
            else:
                comments = ("\n".join(["# {0}".format(c) for c in comments]) +
                            "\n")
        contents = "#%yanny\n" + comments
        #
        # Print any key/value pairs
        #
        for key in self.pairs():
            contents += "{0} {1}\n".format(key, self[key])
        #
        # Print out enum definitions
        #
        if len(self._symbols['enum']) > 0:
            contents += "\n" + "\n\n".join(self._symbols['enum']) + "\n"
        #
        # Print out structure definitions
        #
        if len(self._symbols['struct']) > 0:
            contents += "\n" + "\n\n".join(self._symbols['struct']) + "\n"
        contents += "\n"
        #
        # Print out the data tables
        #
        for sym in self.tables():
            columns = self.columns(sym)
            for k in range(self.size(sym)):
                line = list()
                line.append(sym)
                for col in columns:
                    if self.isarray(sym, col):
                        datum = ('{' + ' '.join([self.protect(x)
                                 for x in self[sym][col][k]]) + '}')
                    else:
                        datum = self.protect(self[sym][col][k])
                    line.append(datum)
                contents += "{0}\n".format(' '.join(line))
        #
        # Actually write the data to file
        #
        with open(newfile, 'w') as f:
            f.write(contents)
        self._contents = contents
        self.filename = newfile
        self._parse()
        return

    def append(self, datatable):
        """Appends data to an existing FTCL/yanny file.

        Tries as much as possible to preserve the ordering & format of the
        original file.  The datatable should adhere to the format of the
        yanny object, but it is not necessary to reproduce the 'symbols'
        dictionary.  It will not try to append data to a file that does not
        exist.  If the append is successful, the data in the object will be
        updated.

        Parameters
        ----------
        datatable : :class:`dict`
            The data to append.
        """
        if len(self.filename) == 0:
            raise ValueError("No filename is set for this object. " +
                             "Use the filename attribute to set the filename!")
        if not isinstance(datatable, dict):
            raise ValueError("Data to append is not of the correct type. " +
                             "Use a dict!")
        timestamp = datetime.datetime.utcnow().strftime(
                    '%Y-%m-%d %H:%M:%S UTC')
        contents = ''
        #
        # Print any key/value pairs
        #
        for key in datatable.keys():
            if key.upper() in self.tables() or key == 'symbols':
                continue
            contents += "{0} {1}\n".format(key, datatable[key])
        #
        # Print out the data tables
        #
        for sym in self.tables():
            if sym.lower() in datatable:
                datasym = sym.lower()
            else:
                datasym = sym
            if datasym in datatable:
                columns = self.columns(sym)
                for k in range(len(datatable[datasym][columns[0]])):
                    line = list()
                    line.append(sym)
                    for col in columns:
                        if self.isarray(sym, col):
                            datum = ('{' + ' '.join([self.protect(x)
                                     for x in datatable[datasym][col][k]]) +
                                     '}')
                        else:
                            datum = self.protect(datatable[datasym][col][k])
                        line.append(datum)
                    contents += "{0}\n".format(' '.join(line))
        #
        # Actually write the data to file
        #
        if len(contents) > 0:
            contents = ("# Appended by yanny.py at {0}.\n".format(timestamp) +
                        contents)
            if os.access(self.filename, os.W_OK):
                with open(self.filename, 'a') as f:
                    f.write(contents)
                self._contents += contents
                self._parse()
            else:
                raise PydlutilsException(self.filename +
                                         " does not exist, aborting append!")
        else:
            warnings.warn("Nothing to be appended!", PydlutilsUserWarning)
        return

    def _parse(self):
        """Converts text into tables that users can use.

        This method is for use internally by the yanny object.  It is not
        meant to be called by users.

        Parsing proceeds in this order:

        #. Lines that end with a backslash character ``\`` are reattached
           to following lines.
        #. Structure & enum definitions are identified, saved into the
           'symbols' dictionary & stripped from the contents.
        #. Structure definitions are interpreted.
        #. At this point, the remaining lines of the original file can only
           contain these things:

           * 'blank' lines, including lines that only contain comments
           * keyword/value pairs
           * structure rows

        #. The remaining lines are scanned sequentially.

           #. 'Blank' lines are identified & ignored.
           #. Whitespace & comments are stripped from non-blank lines.
           #. Empty double braces ``{{}}`` are converted into empty double
              quotes ``""``.
           #. If the first word on a line matches the name of a structure,
              the line is broken up into tokens & each token or set of tokens
              (for arrays) is converted to the appropriate Python type.
           #. If the first word on a line does not match the name of a
              structure, it must be a keyword, so this line is interpreted
              as a keyword/value pair.  No further processing is done to
              the value.

        #. At the conclusion of parsing, if ``self.raw`` is ``False``, the
           structures are converted into NumPy record arrays.
        """
        #
        # there are five things we might find
        # 1. 'blank' lines including comments
        # 2. keyword/value pairs (which may have trailing comments)
        # 3. enumeration definitions
        # 4. structure definitions
        # 5. data
        #
        lines = self._contents
        #
        # Reattach lines ending with \
        #
        lines = re.sub(r'\\\s*\n', ' ', lines)
        #
        # Find structure & enumeration definitions & strip them out
        #
        self._symbols['struct'] = re.findall(
                                    r'typedef\s+struct\s*\{[^}]+\}\s*\w+\s*;',
                                    lines)
        self._symbols['enum'] = re.findall(
                                  r'typedef\s+enum\s*\{[^}]+\}\s*\w+\s*;',
                                  lines)
        lines = re.sub(r'typedef\s+struct\s*\{[^}]+\}\s*\w+\s*;', '', lines)
        lines = re.sub(r'typedef\s+enum\s*\{[^}]+\}\s*\w+\s*;', '', lines)
        #
        # Interpret the structure definitions
        #
        typedefre = re.compile(r'typedef\s+struct\s*\{([^}]+)\}\s*(\w*)\s*;')
        for typedef in self._symbols['struct']:
            typedefm = typedefre.search(typedef)
            (definition, name) = typedefm.groups()
            self[name.upper()] = dict()
            self._symbols[name.upper()] = list()
            definitions = re.findall(r'\S+\s+\S+;', definition)
            for d in definitions:
                d = d.replace(';', '')
                (datatype, column) = re.split(r'\s+', d)
                column = re.sub(r'[\[<].*[\]>]$', '', column)
                self._symbols[name.upper()].append(column)
                self[name.upper()][column] = list()
        # Remove lines containing only comments
        comments = re.compile(r'^\s*#')
        # Remove lines containing only whitespace
        blanks = re.compile(r'^\s*$')
        #
        # Remove trailing comments, but not if they are enclosed in quotes.
        #
        # trailing_comments = re.compile(r'\s*\#.*$')
        # trailing_comments = re.compile(r'\s*\#[^"]+$')
        #
        # Double empty braces get replaced with empty quotes
        #
        double_braces = re.compile(r'\{\s*\{\s*\}\s*\}')
        if len(lines) > 0:
            for line in lines.split('\n'):
                if len(line) == 0:
                    continue
                if comments.search(line) is not None:
                    continue
                if blanks.search(line) is not None:
                    continue
                #
                # Remove leading & trailing blanks & comments
                #
                line = line.strip()
                line = self.trailing_comment(line)
                # line = trailing_comments.sub('',line)
                line = double_braces.sub('""', line)
                #
                # Now if the first word on the line does not match a
                # structure definition it is a keyword/value pair
                #
                (key, value) = self.get_token(line)
                uckey = key.upper()
                if uckey in self._symbols:
                    #
                    # Structure data
                    #
                    for column in self._symbols[uckey]:
                        if len(value) > 0 and blanks.search(value) is None:
                            (data, value) = self.get_token(value)
                            if self.isarray(uckey, column):
                                #
                                # An array value
                                # if it's character data, it won't be
                                # delimited by {} unless it is a multidimensional
                                # string array.  It may or may not be delimited
                                # by double quotes
                                #
                                # Note, we're assuming here that the only
                                # multidimensional arrays are string arrays
                                #
                                arraydata = list()
                                while len(data) > 0:
                                    (token, data) = self.get_token(data)
                                    arraydata.append(token)
                                self[uckey][column].append(
                                    self.convert(uckey, column, arraydata))
                            else:
                                #
                                # A single value
                                #
                                self[uckey][column].append(
                                    self.convert(uckey, column, data))
                        else:
                            break
                else:
                    #
                    # Keyword/value pair
                    #
                    self[key] = value
        #
        # If self.raw is False, convert tables into NumPy record arrays
        #
        if not self.raw:
            for t in self.tables():
                record = np.zeros((self.size(t),), dtype=self.dtype(t))
                for c in self.columns(t):
                    record[c] = self[t][c]
                self[t] = record.view(np.recarray)
        return


def write_ndarray_to_yanny(filename, datatables, structnames=None,
                           enums=None, hdr=None, comments=None):
    """Converts a NumPy record array into a new FTCL/yanny file.

    Returns a new yanny object corresponding to the file.

    Parameters
    ----------
    filename : :class:`str`
        The name of a parameter file.
    datatables : :class:`numpy.ndarray`, :class:`numpy.recarray` or :class:`list` of these.
        A NumPy record array containing data that can be copied into a
        `yanny` object.
    structnames : :class:`str` or :class:`list` of :class:`str`, optional
        The name(s) to give the structure(s) in the yanny file.  Defaults to
        'MYSTRUCT0'.
    enums : :class:`dict`, optional
        A dictionary containing enum information.  See the documentation for
        the :meth:`~pydl.pydlutils.yanny.yanny.dtype_to_struct` method of the
        yanny object.
    hdr : :class:`dict`, optional
        A dictionary containing keyword/value pairs for the 'header' of the
        yanny file.
    comments : :class:`str` or :class:`list` of :class:`str`, optional
        A string containing comments that will be added to the start of the
        new file.

    Returns
    -------
    :class:`yanny`
        The `yanny` object resulting from writing the file.

    Raises
    ------
    PydlutilsException
        If `filename` already exists, or if the metadata are incorrect.
    """
    par = yanny(filename)
    if par:
        #
        # If the file already exists
        #
        raise PydlutilsException(
              "Apparently {0} already exists.".format(filename))
    if isinstance(datatables, (np.ndarray, np.recarray, Table)):
        datatables = (datatables,)
    if structnames is None:
        structnames = ["MYSTRUCT{0:d}".format(k)
                       for k in range(len(datatables))]
    if isinstance(structnames, (str,)):
        structnames = (structnames,)
    if len(datatables) != len(structnames):
        raise PydlutilsException(
              "The data tables and their names do not match!")
    for k in range(len(datatables)):
        struct = par.dtype_to_struct(datatables[k].dtype,
                                     structname=structnames[k], enums=enums)
        par._symbols['struct'] += struct['struct']
        par._symbols[structnames[k].upper()] = struct[structnames[k].upper()]
        if enums is not None and len(par._symbols['enum']) == 0:
            par._symbols['enum'] = struct['enum']
        par[structnames[k].upper()] = datatables[k]
    if hdr is not None:
        for key in hdr:
            par[key] = hdr[key]
    par.write(filename, comments=comments)
    return par


def is_yanny(origin, path, fileobj, *args, **kwargs):
    """Identifies Yanny files or objects.

    This function is for use with
    :func:`~astropy.io.registry.register_identifier`.

    Parameters
    ----------
    origin : :class:`str`
        'read' or 'write'
    path : :class:`str`
        Path to the file.
    fileobj : file object
        Open file object, if available.

    Returns
    -------
    :class:`bool`
        ``True`` if the file or object is a Yanny file.
    """
    if fileobj is not None:
        loc = fileobj.tell()
        fileobj.seek(0)
        try:
            signature = fileobj.read(7)
        finally:
            fileobj.seek(loc)
        return signature == b'#%yanny'
    elif path is not None:
        return path.endswith('.par')
    return isinstance(args[0], yanny)


def read_table_yanny(filename, tablename=None):
    """Read a yanny file into a :class:`~astropy.table.Table`.

    Because yanny files can contain multiple tables, it is necessary
    to specify the table name to return.  However, all "headers"
    (keyword-value pairs) will be included in the Table metadata.

    This function is for use with
    :func:`~astropy.io.registry.register_reader`.

    Parameters
    ----------
    filename : :class:`str`
        Name of the file to read.
    tablename : :class:`str`
        The name of the table to read from the file.

    Returns
    -------
    :class:`~astropy.table.Table`
        The table read from the file.

    Raises
    ------
    :exc:`PydlutilsException`
        If `tablename` is not set.
    :exc:`KeyError`
        If `tablename` does not exist in the file.
    """
    if tablename is None:
        raise PydlutilsException("The name of the table is required!")
    #
    # When opened by Table.read(), the filename will actually be a file-like
    # object opened in *binary* mode.
    #
    par = yanny(filename)
    try:
        t0 = par[tablename.upper()]
    except KeyError:
        raise KeyError("No table named {0} in {1}!".format(tablename, filename))
    t = Table(t0)
    t.meta = par.new_dict_from_pairs()
    return t


def write_table_yanny(table, filename, tablename=None, overwrite=False):
    """Write a :class:`~astropy.table.Table` to a yanny file.

    If `table` has any metadata, it will be written to the file as well.

    This function is for use with
    :func:`~astropy.io.registry.register_writer`.

    Parameters
    ----------
    table : :class:`astropy.table.Table`
        The object to be written.
    filename : :class:`str`
        Name of the file to write to.
    tablename : :class:`str`, optional
        Name to give `table` within the file.
    overwrite : :class:`bool`, optional
        If ``True``, any existing file will be silently overwritten.
    """
    if overwrite and os.path.exists(filename):
        os.remove(filename)
    if table.meta:
        hdr = table.meta
    else:
        hdr = None
    write_ndarray_to_yanny(filename, table, structnames=tablename,
                           hdr=hdr, comments='Table')
