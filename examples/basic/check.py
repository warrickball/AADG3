#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
import AADG3
from astropy.stats import LombScargle
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-N', type=int, default=100,
                    help="number of power spectra to average (default=100)")
args = parser.parse_args()

def LS(t, y):
    f, p = LombScargle(t, y).autopower(normalization='psd',
                                       nyquist_factor=1,
                                       samples_per_peak=1)
    f = f*1e6  # to uHz
    p = p*np.mean(y**2)/np.sum(p)/(f[1]-f[0])  # Bill's normalization
    return f, p

nml, modes, rot = AADG3.load_all_input('basic.in')

y0 = np.loadtxt(nml['output_filename'])
t0 = np.arange(len(y0), dtype=float)*nml['cadence']
y = y0.reshape((args.N, -1))
t = t0[:len(y[0])]

# f, p_tot = LS(t0, y0)
# pl.loglog(f, p_tot)

f, p_tot = LS(t, y[0])

for yi in y[1:]:
    p_tot += LS(t, yi)[1]

pl.loglog(f, p_tot/args.N)
ff = np.linspace(f[0], f[-1], 10000)
pl.loglog(ff, AADG3.PS_model(ff, nml, modes, rot))
pl.show()

# f0, p0 = LS(t0, y0)
# pl.loglog(f0, p0/AADG3.PS_model(f0, nml, modes, rot))
# pl.show()
