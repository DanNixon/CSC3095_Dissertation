#!/usr/bin/env python

# Converts from sigma parameters to kinetic energy for a given mass

import math
import sys
import scipy.constants as scp

def eka(mass, sigma):
    return (scp.value('Planck constant in eV s')**2 * sigma**2) / \
           (2 * mass * scp.value('atomic mass constant') * scp.milli)

if __name__ == '__main__':
    ek = eka(float(sys.argv[1]), float(sys.argv[2]))
    print '{0} meV'.format(ek)
