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
parser.add_argument('-o', '--output-files', type=str, nargs='+', default=None,
                    help="use this timeseries file instead of the one in "
                    "input file")
parser.add_argument('--legend', action='store_true',
                    help="add a legend, using filenames as keys")
args = parser.parse_args()

def LS(t, y):
    f, p = LombScargle(t, y).autopower(normalization='psd',
                                       nyquist_factor=1,
                                       samples_per_peak=1)
    f = f*1e6  # to uHz
    p = p*np.var(y)/np.trapz(p, x=f)  # Bill's normalization
    return f, p

nml, modes, rot = AADG3.load_all_input(args.filename)

filenames = args.output_files if args.output_files else [nml['output_filename']]

for filename in filenames:
    try:
        y0 = np.loadtxt(filename)
    except IOError:
        y0 = np.loadtxt('/'.join(args.filename.split('/')[:-1] + filename))

    t0 = np.arange(len(y0), dtype=float)*nml['cadence']
    y = y0[:len(y0)//args.N*args.N].reshape((args.N, -1))
    t = t0[:len(y[0])]

    f, p_tot = LS(t, y[0])

    for yi in y[1:]:
        p_tot += LS(t, yi)[1]

    pl.loglog(f, p_tot/args.N, label=filename)
    
ff = np.linspace(f[0], f[-1], 10000)
pl.loglog(ff, AADG3.PS_model(ff, nml, modes, rot), label=args.filename)

if args.legend:
    pl.legend()
    
pl.show()

# f0, p0 = LS(t0, y0)
# pl.loglog(f0, p0/AADG3.PS_model(f0, nml, modes, rot))
# pl.show()
