"""Python tools for working with the asteroFLAG Artificial Dataset
Generator."""

import numpy as np

def PS_model(f, nml, modes, rot):
    """Create a mean power spectrum model that roughly matches the AADG3
    input given by namelist ``nml``, mode parameters ``modes`` and
    rotational splittings ``rot`` (requires SciPy).

    The model currently assumes that the power spectra of the modes
    are Lorentzians, therefore ignoring any mode asymmetry introduced
    by the ``rho`` parameter.

    Parameters
    ----------
    f: array
        Frequency range for the power spectrum model
    nml: dict
        Dictionary describing Fortran namelist
    modes: recarray
        Structured array with mode frequency information from ``modes``
        file, as returned by ``load_modes``.
    rot: recarray
        Structured array with rotational splitting information from
        ``rot`` file, as returned by ``load_rot``.

    Returns
    -------
    p: array
        Mean power spectrum at the frequencies given by ``f``.

    """
    from math import factorial
    from scipy.special import lpmn

    x = np.cos(np.radians(nml['inclination']))
    lmax = max(modes['l'])
    E = np.array(lpmn(lmax, lmax, x)[0].T)
    for l in range(lmax+1):
        for m in range(l+1):
            E[l][m] = E[l][m]**2*factorial(l-m)/factorial(l+m)
    
    # De Ridder et al., MNRAS 365, 595 (2006), eq. (19)
    # /2 because AADG works in rms amplitude and /1e6 for ppm^2/uHz
    # note the sinc function for apodisation (np.sinc(x) = sin(pi*x)/(pi*x))
    # and reflection of super-Nyquist part of background
    f6 = f/1e6
    bkg = 1/(1 + (f6*nml['tau']*2.*np.pi)**2)*np.sinc(nml['cadence']*f6)**2
    f6 = (f[::-1] + f[-1])/1e6
    bkg += 1/(1 + (f6*nml['tau']*2.*np.pi)**2)*np.sinc(nml['cadence']*f6)**2
    bkg *= 4.*nml['sig']**2*nml['tau']/2e6

    def lorentz(nu, nu0, w):
        return 1./(1.+4.*(nu-nu0)**2/w**2)

    osc = np.zeros_like(f)
    for row in modes:
        l = row['l']
        height = 2.*row['power']/np.pi/row['width']*nml['p(%i)' % l]
        osc += height*E[l][0]*lorentz(f, row['freq'], row['width'])

        for m in range(1, l+1):
            splitting = rot[(rot['l']==l)
                            &(rot['m']==m)
                            &(rot['n']==row['n'])]['splitting']
            if len(splitting) > 1:
                splitting = splitting[0]

            osc += height*E[l][m]* \
                   (lorentz(f, row['freq']+m*splitting, row['width']) +
                    lorentz(f, row['freq']-m*splitting, row['width']))

    if nml['add_granulation']:
        return bkg + osc
    else:
        return osc
    

def load_all_input(filename):
    """Load all data specified by namelist ``filename``, including mode
    data and rotational splittings in ``modes_filename`` and
    ``rotation_filename``.  If ``rotation_filename`` is ``''``, return
    rotation data equivalent to zero rotation (which is what the
    simulator does).

    Parameters
    ----------
    filename: str
        Filename of the main namelist (``&controls``).

    Returns
    -------
    nml: dict
        Dictionary describing Fortran namelist
    modes: recarray
        Structured array with mode frequency information from ``modes``
        file, as returned by ``load_modes``.
    rot: recarray
        Structured array with rotational splitting information from
        ``rot`` file, as returned by ``load_rot``.

    """
    nml = load_namelist(filename)
    base = '/'.join(filename.split('/')[:-1]) + '/'
    try:
        modes = load_modes(nml['modes_filename'])
    except OSError:
        modes = load_modes(base + nml['modes_filename'])

    if nml['rotation_filename'] == '':
        rot = generate_const_rot(modes)
    else:
        try:
            rot = load_rot(nml['rotation_filename'])
        except OSError:
            rot = load_rot(base + nml['rotation_filename'])

    return nml, modes, rot


def load_namelist(filename):
    """Load input data from ``controls`` namelist in ``filename``.
    Returns a Python dict in which each namelist variable is a key
    with the appropriate value. e.g. the namelist
    ::

        &controls
            n_fine = 50
        /

    becomes the Python dict
    ::

        {'n_fine': 50}

    """
    # include valid default values
    nml = {'user_seed': 0, 'n_fine': 0, 'inclination': 0.0,
           'cycle_period': 1e99, 'cycle_phase': 0.0, 'nuac': 0.0, 'sdnu': 0.0,
           'p(0)': 1.0, 'p(1)': 0.0, 'p(2)': 0.0, 'p(3)': 0.0,
           'add_granulation': True, 'rotation_filename': '',
           'verbose': False}

    with open(filename, 'r') as f:
        lines = [[word.strip() for word in line.split('=')[:2]] for
                 line in f.readlines() if '=' in line and not
                 line.strip().startswith('#') and not
                 line.strip().startswith('!')]

    # strip trailing comments
    for line in lines:
        try:
            line[1] = line[1][:line[1].index('!')].strip()
        except ValueError:
            continue

    # parse input data into dict
    for line in lines:
        k = line[0].lower()
        v = line[1]
        if v == '.true.':
            nml[k] = True
        elif v == '.false.':
            nml[k] = False
        elif v[0] == "'" and v[-1] == "'":
            nml[k] = v[1:-1]
        elif '.' in v or 'd' in v.lower() or 'e' in v.lower():
            nml[k] = float(v.replace('d','e'))
        else:
            nml[k] = int(v)

    return nml


modes_dtype = [('l', int), ('n', int), ('freq', float),
               ('width', float), ('power', float), ('dfreq', float)]

def load_modes(filename):
    """Load data for oscillation modes from ``filename``.  Returns a NumPy
    record array with fields corresponding to:

    * angular degree (``l``),
    * radial order (``n``),
    * mode frequency in uHz (``freq``),
    * mode linewidth in uHz (``width``),
    * mode power (``power``) and
    * activity cycle variation (``dfreq``).
    
    """
    return np.loadtxt(filename, dtype=modes_dtype, comments=['#', '!'])

rot_dtype = [('n', int), ('l', int), ('m', int),
             ('splitting', float)]

def load_rot(filename):
    """Load data for rotational splittings from ``filename``.  Returns a
    NumPy record array with fields corresponding to:

    * radial order (``n``),
    * angular degree (``l``),
    * azimuthal order (``m``) and
    * the rotational splitting (``splitting``).

    Note that the AADG3's input rotational splitting is ``1/m`` times
    the actual frequency difference from ``m=0`` in the power
    spectrum. i.e. the mode frequency of an ``m=2`` mode is the radial
    mode frequency plus twice the splitting in the output array.

    """
    return np.loadtxt(filename, dtype=rot_dtype, comments=['#', '!'])


def generate_const_rot(modes, splitting=0.0):
    """Create correctly-formatted data for rotational splittings with a
    constant value (which is by default zero).  Useful if you have
    mode data but not rotation data.

    Parameters
    ----------
    modes: recarray
        Structured array with mode frequency information from ``modes``
        file, as returned by ``load_modes``.
    splitting: float, optional
        Constant rotational splitting in uHz. (default=0)
    
    Returns
    -------
    rot: recarray
        Structured array with rotational splitting information in the
        same format as returned by ``load_rot``.

    """
    rot = np.zeros(sum(modes['l']), dtype=rot_dtype)
    rot['splitting'] = splitting
    i = 0
    for row in modes:
        for m in range(1, row['l']+1):
            rot['n'][i] = row['n']
            rot['l'][i] = row['l']
            rot['m'][i] = m
            i += 1
            
    return rot
    

def save_namelist(filename, nml):
    """Saves the Fortran namelist in Python dict ``nml`` to a file named
    ``filename``.  The Python dict should have keys corresponding to
    the namelist variable names, as returned by ``load_nml``. e.g. the
    Python dict
    ::

        {'n_fine': 50}

    would be saved as
    ::

        &controls
            n_fine = 50
        /
    """
    with open(filename, 'w') as f:
        f.write('&controls\n')
        for k, v in nml.items():
            if type(v) == bool:
                if v:
                    f.write('    %s = .true.\n' % k)
                else:
                    f.write('    %s = .false.\n' % k)
            elif type(v) == float or type(v) == np.float or type(v) == np.float64:
                f.write('    %s = %gd0\n' % (k, v))
            elif type(v) == int:
                f.write('    %s = %i\n' % (k, v))
            else:
                f.write("    %s = '%s'\n" % (k, v))

        f.write('/\n')

def new_namelist(mission=None):
    nml = {'p(0)': 1.0, 'p(1)': 1.5, 'p(2)': 0.5, 'p(3)': 0.05}

    # there are 1461/4 days per year
    if mission.lower() == 'tess':
        # 27.4d at 2min cadence
        nml['cadence'] = 120.0
        nml['n_cadences'] = 720*137//5
    elif mission.lower() == 'kepler_sc':
        # 4yr at 1min cadence
        nml['cadence'] = 60.0
        nml['n_cadences'] = 4*1440*1461//4  # int(4.0*365.25*1440.0)
    elif mission.lower() == 'kepler_lc':
        # 4yr at 30min cadence
        nml['cadence'] = 60.0
        nml['n_cadences'] = 4*24*2*1461//4 # int(4.0*365.25*24.0*2.0)
    elif mission.lower() == 'k2_sc':
        # 90d at 1min cadence
        nml['cadence'] = 60.0
        nml['n_cadences'] == 90*1440
    elif mission.lower() == 'k2_lc':
        # 90d at 30min cadence
        nml['cadence'] = 60.0
        nml['n_cadences'] == 90*24*2
        
    return nml
