#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
import AADG3
from astropy.stats import LombScargle
from argparse import ArgumentParser

parser = ArgumentParser(description="""
Reads an AADG3 input file and compares a (possibly averaged) power
spectrum to the expected (though simplified) model power spectrum.
The model power spectrum currently doesn't include mode asymmetry.
""")
parser.add_argument('filename', type=str,
                    help="AADG3 input filename")
parser.add_argument('-N', type=int, default=1,
                    help="number of power spectra to average (default=1, "
                    "i.e. no averaging)")
args = parser.parse_args()

def LS(t, y):
    f, p = LombScargle(t, y).autopower(normalization='psd',
                                       nyquist_factor=1,
                                       samples_per_peak=1)
    f = f*1e6  # to uHz
    p = p*np.mean(y**2)/np.sum(p)/(f[1]-f[0])  # Bill's normalization
    return f, p

nml, modes, rot = AADG3.load_all_input(args.filename)

y0 = np.loadtxt(nml['output_filename'])
t0 = np.arange(len(y0), dtype=float)*nml['cadence']
y = y0[:len(y0)//args.N*args.N].reshape((args.N, -1))
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
