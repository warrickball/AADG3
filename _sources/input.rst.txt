Input files and command-line arguments
======================================

Namelist control file
+++++++++++++++++++++

- ``user_seed`` (*integer*)

  Seed used for random number generator to facilitate reproducible
  results.

  If greater than zero, the value is used to create a series of
  integers passed to the intrinsic ``random_seed`` function (see the
  code for details).

  If equal to zero, then the random number generator (RNG) is
  initialised by calling ``random_seed`` without an argument.  Note
  that this gives different behaviour for ``gfortran`` major versions
  up to and including 6 versus major versions from 7 onwards.  In the
  former (major versions :math:`\leq6`), the RNG is initialised to
  some default state, which will be the same each time the program is
  run (on a given system).  In the latter (major versions
  :math:`\geq7`), the RNG is initialised using random data from the
  operating system and will produce different output each time the
  program is run.

- ``cadence`` (*float*)

  Cadence of the timeseries in seconds.

- ``n_cadences`` (*integer*)

  Number of cadences in the timeseries.

- ``n_relax`` (*integer*)

  The number of cadences over which the stochastic driving term is
  allowed to build up from zero.

- ``n_fine`` (*integer*)

  When the driving terms ("kicks") are generated, the cadence is
  subdivided into ``n_fine`` subcadences.

- ``sig`` (*float*)

  Rms amplitude of the granulation.  The power spectral density of the
  granulation at frequency :math:`\nu` will be approximately
  :math:`2\sigma^2\tau/(1+(2\pi\nu\tau)^2)`, where :math:`\tau` is the
  granulation timescale as set by ``tau``, above.  A typical solar
  value for ``sig`` is ``60d0``.

- ``rho`` (*float*)

  Relative weighting of the correlated kicks (:math:`u_c`) and
  uncorrelated kicks (:math:`u_u`).  The total kick is :math:`u=\rho
  u_c + \sqrt{1-\rho^2}u_u`, so that ``rho=0`` means purely
  uncorrelated kicks and ``rho=1`` means purely correlated kicks.  See
  `Toutain et al. (2006)
  <http://adsabs.harvard.edu/abs/2006MNRAS.371.1731T>`_ for details.

- ``tau`` (*float*)

  Timescale for the granulation in seconds.  A typical solar value is
  for ``tau`` is about ``250d0``.

- ``inclination`` (*float*)

  Inclination angle of the rotation axis in degrees.  That is, ``90d0``
  is edge-on and ``0d0`` is pole-on.

- ``p(:)`` (*float*, *default p(0)=1.0d0, p(1:)=0.0d0*)

  Visibility ratios for modes of different angular degrees, as
  indicated by the index of the array.  By definition, ``p(0)`` should
  be 1, which is the default.

- ``add_granulation`` (*logical*)

  If ``.true.``, add the granulation signal to the timeseries for the
  oscillation modes.

- ``modes_filename`` (*string*)

  Name of the file containing the mode frequency information (usually
  ending ``.con``).

- ``rotation_filename`` (*string*)

  Name of the file containing the rotational splitting information
  (usually ending ``.rot``).

- ``output_filename`` (*string*)

  Name of the file for the output timeseries.

- ``output_fmt`` (*string*, *default='(f16.7)'*)

  Format statement for the timeseries.

- ``verbose`` (*logical*, *default=.false.*)

  If ``.true.``, print some information to the terminal as the
  simulation progresses.


Mode data file format (``modes_filename``)
++++++++++++++++++++++++++++++++++++++++++

The file ``modes_filename``, specified in the namelist control file,
contains details about the modes that are determined by the radial
order :math:`n` and angular degree :math:`\ell` but not the azimuthal
order :math:`m`. i.e. everything except the rotational splittings.  The file
must have six columns, giving the

- angular degree :math:`\ell`,
- radial order :math:`n`,
- frequency :math:`\nu` in :math:`\mu\mathrm{Hz}`,
- linewidth :math:`\Gamma` in :math:`\mu\mathrm{Hz}`,
- rms power in units of the desired output
  (e.g. :math:`(\mathrm{cm}/\mathrm{s})^2` or :math:`\mathrm{ppm}^2`),
  and
- activity-induced frequency shift in :math:`\mu\mathrm{Hz}`

of each mode.  For example, the line specifying an
:math:`(n,\ell)=(21,1)` mode with no cycle shift, frequency
:math:`3\,\mathrm{mHz}`, linewidth :math:`1\mu\mathrm{Hz}` and power
:math:`10\,\mathrm{ppm}^2` would read something like

::

   1  21   3000.0   1.0   10.0   0.0

Lines starting with ``#`` or ``!`` are ignored.


Rotational splitting file format (``rotation_filename``)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The file ``rotation_filename``, specified in the namelist control
file, contains details about the rotational splittings of the modes,
which depend on the radial order :math:`n`, angular degree
:math:`\ell` and azimuthal order :math:`m`.  The file must have
four columns, giving the

- radial order :math:`n`,
- angular degree :math:`\ell`,
- azimuthal order :math:`m`, and
- rotational splitting :math:`\delta\nu_{nlm}`

of each mode.

The rotation splitting is assumed to be the same for :math:`m` and
:math:`-m`, so the file must only specifiy splittings for positive
:math:`m`.

The rotational splittings will be multiplied by :math:`m`, so an
:math:`m=2` mode with rotational splitting :math:`400\,\mathrm{nHz}`
will be separated from the :math:`m=0` mode by
:math:`800\,\mathrm{nHz}`.  If :math:`(n,\ell)=(21,1)`, this mode
would be specified by the line

::

   21  1  2  0.400

Lines starting with ``#`` or ``!`` are ignored.

  
Command-line arguments
++++++++++++++++++++++

In normal operation, the first command line argument should always be
the namelist control file.  Thereafter, any of the namelist controls
above can be overridden on the command line by invoking AADG3 with

::
   
    AADG3 controls.in --option value

For example, if you wanted to run the timeseries for 40000 cadences
instead of the number the input file, you could use

::
   
    AADG3 controls.in --n_cadences 40000

or

::
   
    AADG3 controls.in --n-cadences 40000

For a boolean (true/false) option, adding it as a command-line
argument sets it to true.  Using ``--no-option`` instead of
``--option`` will set it to false. e.g.

::

   AADG3 controls.in --add-granulation

will set ``add_granulation`` to true, whereas

::

   AADG3 controls.in --no-add-granulation

will set ``add_granulation`` to false.  The ``verbose`` option is a
special case.  It can be set and unset as above but, in addition, it
can be set with ``-v`` and unset with ``--quiet`` or ``-q``.
    
Failure to parse a genuine option is a bug and should be reported.

Other command line options are:

- ``-h``, ``--help``

  Shows some help, then quits.

- ``-V``, ``--version``

  Shows the version number, then quits.

Both of these options must be the first argument after ``AADG3`` if
they are used.
