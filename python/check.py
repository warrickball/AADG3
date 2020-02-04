#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
import AADG3
from astropy.timeseries import LombScargle
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

nml, modes, rot = AADG3.load_all_input(args.filename)

filenames = args.output_files if args.output_files else [nml['output_filename']]

for filename in filenames:
    try:
        y0 = np.loadtxt(filename)
    except IOError:
        y0 = np.loadtxt('/'.join(args.filename.split('/')[:-1] + [filename]))

    y = y0[:len(y0)//args.N*args.N].reshape((args.N, -1))

    f = np.fft.rfftfreq(len(y[0]), d=nml['cadence'])
    f *= 1e6
    p = np.zeros_like(f)

    for yi in y:
        pi = np.abs(np.fft.rfft(yi))**2
        pi = pi*np.var(yi)/np.trapz(pi, x=f) # Parseval's theorem
        p += pi

    pl.loglog(f, p/args.N, label=filename)

# create a frequency range that resolves the analytic Lorentzians
ff = []
one = np.tan(np.linspace(-np.pi/2.2, np.pi/2.2, 21))
one = one/np.max(one)*5
for row in modes:
    l = row['l']

    ff.append(one*row['width'] + row['freq'])

    for m in range(1, l+1):
        splitting = rot[(rot['l']==l)
                        &(rot['m']==m)
                        &(rot['n']==row['n'])]['splitting']
        if len(splitting) > 1:
            splitting = splitting[0]

        ff.append(one*row['width'] + row['freq'] - m*splitting)
        ff.append(one*row['width'] + row['freq'] + m*splitting)

# combine frequency mesh for modes with 5000 points across background
ff = np.sort(np.hstack(ff + [np.linspace(f[0], f[-1], 5000)]))

pl.loglog(ff, AADG3.PS_model(ff, nml, modes, rot), label=args.filename)

if args.legend:
    pl.legend()

pl.show()

# f0, p0 = LS(t0, y0)
# pl.loglog(f0, p0/AADG3.PS_model(f0, nml, modes, rot))
# pl.show()
