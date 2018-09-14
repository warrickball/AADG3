.. AADG3 documentation master file, created by
   sphinx-quickstart on Fri Mar 23 17:01:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The AsteroFLAG Artificial Dataset Generator 3 (AADG3)
=====================================================

Introduction
++++++++++++

AADG3 simulates lightcurves of solar-like oscillators.  It is based on
a series of programs derived from the original work by `Chaplin et
al. (1997)`_, which describes the method of generating the
oscillations using the Laplace transform solution of a damped harmonic
oscillator.  The use of a correlated excitation component to produce
mode asymmetry was introduced by `Toutain et al. (2006)`_.  Most
recently, `Howe et al. (2015)`_ gave a detailed description of the
code's methods as part of their validation of BiSON results.

Variants of the code have been used by the *Solar Fitting at Low
Angular Degree* (solarFLAG, e.g. `Chaplin et al. 2006`_) and
subsequent *Asteroseismic Fitting at Low Angular Degree* (asteroFLAG,
e.g. `Chaplin et al. 2008`_) consortia.

This version is derived from code delivered as version 2, which we
incremented to the present version.  We publicly released the code as
part of Ball et al. (in prep.), in which we also presented a large
catalogue of mock TESS targets for which we generated data using
AADG3.

Please cite these papers appropriately if you use AADG3 in your
research.

.. _`Chaplin et al. (1997)`: http://adsabs.harvard.edu/abs/1997MNRAS.287...51C
.. _`Chaplin et al. 2006`: http://adsabs.harvard.edu/abs/2006MNRAS.369..985C
.. _`Toutain et al. (2006)`: http://adsabs.harvard.edu/abs/2006MNRAS.371.1731T
.. _`Chaplin et al. 2008`: http://adsabs.harvard.edu/abs/2008AN....329..549C
.. _`Howe et al. (2015)`: http://adsabs.harvard.edu/abs/2015MNRAS.454.4120H

Download
++++++++

AADG3 is hosted on `GitHub <https://github.com/warrickball/AADG3>`_, where you can find and download from the 
`list of release versions <https://github.com/warrickball/AADG3/releases>`_ or fork/clone
the repo.  For convenience, from here you can download archives of the

* most recent release, **v3.0.0a** (`tarball <https://github.com/warrickball/AADG3/archive/v3.0.0a.tar.gz>`__,
  `zip <https://github.com/warrickball/AADG3/archive/v3.0.0a.zip>`__),
  
* current stable branch (`tarball <https://github.com/warrickball/AADG3/archive/stable.tar.gz>`__, 
  `zip <https://github.com/warrickball/AADG3/archive/stable.zip>`__), and
  
* current development version (master branch, `tarball <https://github.com/warrickball/AADG3/archive/master.tar.gz>`__, `zip <https://github.com/warrickball/AADG3/archive/master.zip>`__)
  

Installation
++++++++++++

Extract the archive, change into the directory it creates, then run

::
   
    cd make
    make

This will place the standalone executable at ``bin/AADG3`` under the
AADG3 folder.

The ``makefile`` is set up to compile using ``gfortran``.  Any
reasonably up-to-date version of ``gfortran`` should be fine.  The
code has been partially tested using versions 4.8.5 and 7.2.0 and
more extensively tested with versions 5.4.0 and 8.1.1.  If you're
having problems, you can install an up-to-date version of ``gfortran``
on Mac OS or Linux by using the `MESA SDK
<http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk>`_ or `MAD
SDK <http://www.astro.wisc.edu/~townsend/static.php?ref=madsdk>`_.

Note a subtle but important change of behaviour in the random number
generator (RNG) from ``gfortran`` version 6 to 7.  Before version 7,
calling ``random_seed`` with no input sets the RNG to some default
state, which is always the same (at least on a given system).  From
version 7 onwards, calling ``random_seed`` with no input sets the RNG
using data from the operating system, which is different each time the
program is run.  Thus, if you set ``user_seed=0`` when compiling with
versions before 7, you'll get the same timeseries each time you call
the program.  If you used versions after 7, you'll get a different
timeseries each time.

Running AADG3
+++++++++++++

The standalone executable file is placed in ``bin/AADG3``.  You can
run it from anywhere on your system or add it to your ``$PATH``
environment variable.  I personally add a symbolic link to my
``~/.local/bin`` folder, which is in my ``$PATH``.

You can run some basic tests using ``bin/test_AADG3``.  Please let me
know if any failures are reported.

To get started, a few examples are provided in ``examples/``.
Navigate to any of those and run ``AADG3 <example>.in``.  e.g.

::

   cd examples/s4tess_llrgb
   AADG3 llrgb.in

The file ``llrgb.asc`` then contains the output timeseries.  The input
file ``llrgb.in`` contains some information about this example.  These
runs can take several minutes to finish, so don't be worried if you
don't see any output for a little while.

The various inputs are described in detail `here <input.html>`__.

Python module
+++++++++++++

The AADG3 distribution includes a small Python module under
``python/``.  This provides a handful of functions for reading and
writing input and output for the program.  It also has a function
(``PS_model``) to generate an approximate model of the power spectrum
that the code should produce, which is useful for basic testing and
cross-checking of results.  To learn more about the Python module,
please read the `API documentation <python.html>`__.

Reference
+++++++++

.. toctree::
   :maxdepth: 2

   input
   fortran
   python


..
   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
