===========================
ABOUT THE COOLING FUNCTIONS
===========================

Introduction
------------

This file has been adapted from http://www.mso.anu.edu.au/~ralph/data/cool/ABOUT.txt.

The associated files contain the cooling functions calculated for the paper
by `Sutherland & Dopita (1993) <http://adsabs.harvard.edu/abs/1993ApJS...88..253S>`_.
See paper or preprint for details.

FILE LIST
---------

Collisonal Ionisation Equilibrium cooling functions as a function
of metallicity are given in the files:

mzero.cie
    cie cooling for primordial hydrogen/helium mix log(T) = 4-8.5

m-00.cie
    cie cooling for solar abundances mix log(T) = 4-8.5

m-05.cie
    [Fe/H] = -0.5, solar/primordial average ratios

m-10.cie
    [Fe/H] = -1.0, primordial ratios (ie enhanced oxygen)

m-15.cie
    [Fe/H] = -1.5, primordial ratios

m-20.cie
    [Fe/H] = -2.0, primordial ratios

m-30.cie
    [Fe/H] = -3.0, primordial ratios

m+05.cie
    [Fe/H] = +0.5, solar ratios log(T) = 4.1-8.5 (due to charge
    exchange problems at log(T) = 4.0)

FILE FORMAT
-----------

The data files consist of tab separated ASCII column data described as follows:

log(T)
    log temperature in K

ne, nH, nt
    number densities, electrons, hydrogen and total ion in cm^-3.

log(lambda net), log(lambda norm)
    log on the net cooling function and the normalised cooling function.
    lambda norm = lambda net / (ne nt). lambda net in ergs cm^-3 s^-1,
    lambda net in ergs cm^3 s^-1.  While the
    densities are kept less than about p/k 10^8 both isobaric and isochoric
    curves can be constructed from the normalised function using an appropriate
    density function.  The non-equilibrium net cooling function is from the
    isobaric model used to calculate the curves.  In the CIE curves the net function
    is for the isochoric cie model.

log(U)
    U = 3/2 N kT, N = ne + nt the total internal energy. ergs cm^-3

log(taucool)
    The normalized cooling timescale Nr*( U/(lambda net))
    Nr = (ne nt)/(ne+nt).   s cm^-3

P12
    Gas pressure NkT. ergs ergs cm^-3 times 10^12

rho24
    Density g  cm^-3 times 10^24

Ci
    The isothermal sound speed kms s^-1

mubar
    mean molecular weight grams times 10^24

NOTES
-----

Tests showed that density quenching begins to affect the curves in the low
temperature parts of the cooling functions when log(p/k) = 8 or more.

Above log(T) = 7.5-8.5 it should be pretty safe to use a powerlaw
fit to the free-free losses

Good luck, and let me know if you have any problems.

These files will be available via anonymous ftp from

Yours sincerely

Ralph S. Sutherland  1992-1993.
