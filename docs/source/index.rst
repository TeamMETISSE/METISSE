.. METISSE documentation master file, created by
   sphinx-quickstart on Sat May 25 10:38:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to METISSE's documentation!
===================================

**METISSE** is an acronym for **MET**\ hod of **I**\ nterpolation for **S**\ ingle **S**\ tar **E**\ volution. It is a code package designed for rapidly computing the evolution of stellar populations by interpolating within a set of pre-computed evolutionary tracks. The input tracks for interpolation are computed from stellar structure and evolution codes such as MESA. The input tracks must be in EEP format, meaning that significant evolutionary points like zero-age main sequence (ZAMS) should occur at the same line number across each file.  METISSE can be used as a standalone code or in conjunction with other codes for evolving populations of single and binary stars. Its special structure allows it to be easily used in codes that already use the fitting formulae-based code SSE (
`Hurley et al. 2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_ 
).


**Prequisite:** METISSE requires gcc/6.4.0 or above

The code package is available `here <https://github.com/TeamMETISSE/METISSE>`_.

In addition to the code, METISSE also requires a set of EEP tracks with a given metallicity to interpolate a stellar track for the same metallicity and a specific mass.
EEP tracks for a range of mass and metallicity, computed using MESA, and ready for use with METISSE can be downloaded from this link. 

After downloading, the path of this folder should be provided to METISSE. This path and other inputs should be provided to METISSE depending on how METISSE is used (in standalone mode or in conjunction with other codes).

More details about METISSE including code capabilities are described in:

- `Agrawal, Poojan; Hurley, Jarrod; Stevenson, Simon; Szécsi, Dorottya; Flynn, Chris 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.4549A/abstract>`_
- `Agrawal, Poojan; Hurley, Jarrod; Stevenson; Rodriguez, Carl L.; Szécsi, Dorottya; Kemp, Alex 2023 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..933A/abstract>`_

Check out the following sections for further information.


.. toctree::
   :maxdepth: 0

   contents
   usage_standalone
   usage_other
   eep_tracks

