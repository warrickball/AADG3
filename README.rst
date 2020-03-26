The AsteroFLAG Artificial Data Generator 3 (AADG3)
==================================================

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
incremented to the present version.  We publicly released version 3.0 as
part of `Ball et al. (2018)`_, in which we also presented a large
catalogue of mock TESS targets for which we generated data using
AADG3.

Please cite these papers appropriately if you use AADG3 in your
research.

.. _`Chaplin et al. (1997)`: https://ui.adsabs.harvard.edu/abs/1997MNRAS.287...51C/abstract
.. _`Chaplin et al. 2006`: https://ui.adsabs.harvard.edu/abs/2006MNRAS.369..985C/abstract
.. _`Toutain et al. (2006)`: https://ui.adsabs.harvard.edu/abs/2006MNRAS.371.1731T/abstract
.. _`Chaplin et al. 2008`: https://ui.adsabs.harvard.edu/abs/2008AN....329..549C/abstract
.. _`Howe et al. (2015)`: https://ui.adsabs.harvard.edu/abs/2015MNRAS.454.4120H/abstract
.. _`Ball et al. (2018)`: https://ui.adsabs.harvard.edu/abs/2018ApJS..239...34B/abstract

Documentation
-------------

Documentation for AADG3 is available `here <https://warrickball.github.io/AADG3/>`__.


License
-------

Copyright 2018-2020 Warrick Ball & Bill Chaplin

AADG3 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AADG3 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AADG3.  If not, see `<https://www.gnu.org/licenses/>`_.
