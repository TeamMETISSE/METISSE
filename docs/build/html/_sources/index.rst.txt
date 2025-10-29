.. METISSE documentation master file, created by
   sphinx-quickstart on Sat May 25 10:38:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to METISSE's documentation!
====================================


**METISSE** is an acronym for **MET**\ hod of **I**\ nterpolation for **S**\ ingle **S**\ tar **E**\ volution. It is an open-source stellar evolution tool designed to provide a fast and flexible alternative to traditional rapid stellar evolution prescriptions.

At its core, METISSE interpolates between pre-computed stellar evolution tracks to determine stellar parameters such as luminosity, radius, and core mass at any evolutionary stage. This approach allows users to incorporate stellar models directly into large-scale stellar and binary population synthesis simulations.


What METISSE is for
-------------------

METISSE has been specifically designed to be integrated with binary evolution and population synthesis frameworks that currently rely on the Fortran 77 based SSE `(Hurley et al. 2000) <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_ codes. METISSE enables such frameworks to:

- Incorporate the latest stellar models from codes like MESA or the BonnCode,
- Vary stellar physics assumptions systematically,
- Explore how uncertainties in stellar evolution impact populations of single and binary stars.

Although METISSE can be run as a standalone program to evolve individual stars or simple stellar populations, its main strength lies in serving as an stellar engine within population synthesis environments.

What you can do with METISSE
----------------------------

- Run METISSE independently to evolve single stellar systems from pre-computed tracks
- Use METISSE with supported binary evolution codes to evolve stars in binaries 
- Use your own custom stellar tracks to test different physical assumptions
- Combine METISSE with population synthesis tools to perform large-scale stellar evolution studies using modern stellar data


.. Learn more
.. ----------

.. - **Installation and setup** — how to download, compile, and prepare METISSE
.. - **Input stellar tracks** — how to use your own models
.. - **Examples** — example runs for single and binary stars
.. - **API and structure** — how METISSE works under the hood

.. Learn more
.. ----------

.. raw:: html

   <div>
   <h3>Learn more</h3>

   <ul>
     <li><a href="installation.html"><b>Installation and setup</b> — how to download METISSE</a></li>
     <li><a href="usage_standalone.html"><b>Single star evolution</b> with METISSE</a></li>
     <li><a href="usage_other.html"><b>Binary evolution and population synthesis</b> with METISSE as the stellar engine</a></li>
     <li><a href="using_custom_input_tracks.html"><b>Input stellar tracks</b> — how to use your own models</a></li>
     <li><a href="structure.html"><b>Structure and Integration</b> — how METISSE works under the hood</a></li>
   </ul>
   </div>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Contents

   installation
   usage_standalone
   usage_other
   input_files
   using_custom_input_tracks
   structure
   FAQ
   glossary

   


.. **METISSE** is an acronym for **MET**\ hod of **I**\ nterpolation for **S**\ ingle **S**\ tar **E**\ volution. It is a code package designed for rapidly computing the evolution of stellar populations by interpolating within a set of pre-computed stellar evolutionary tracks. These input tracks can be computed from various stellar structure and evolution codes.
 
.. It can be used as a `standalone code <usage_standalone.html>`_ for evolving populations of single stars or in conjunction with `binary population synthesis codes <usage_other.html>`_ for evolving populations of single and binary stars. 
.. Check out the following sections for more information about METISSE and how to use it.


..    installation
..    usage_standalone
..    usage_other
..    input_files
..    using_custom_input_tracks
..    structure
..    FAQ
..    glossary
