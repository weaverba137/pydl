#
#
#
def current_mjd():
    import time
    mjd = time.time()/86400.0 + 40587.0
    return mjd

if __name__ == '__main__':
    print(current_mjd())

