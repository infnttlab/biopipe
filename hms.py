import sys

seconds = int(sys.argv[1])

def hms(x):
    m, s = divmod(x, 60)
    h, m = divmod(m, 60)
    return r"%d:%02d:%02d" % (h, m, s)

print (hms(seconds))