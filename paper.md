---
title: 'METISSE: METhod of Interpolation for Single Star Evolution'
tags:
  - Fortran
  - astronomy
  - stellar populations
  - stellar evolution
  - stars
authors:
  - name: Poojan Agrawal
    orcid: 0000-0002-1135-984X
    affiliation: "1,2" # (Multiple affiliations must be quoted)
  - name: Katie Brievik
    orcid: 0000-0001-5228-6598
    affiliation: 3 
  - name: Jarrod Hurley
    orcid: 0000-0003-2694-0415
    affiliation: "4,5" 
  - name: Duncan Maclean
    orcid: 
    affiliation: 2 
  - name: Carl Rodriguez
    orcid: 0000-0003-4175-8881
    affiliation: "2,3" 
  - name: Alex Kemp
    orcid: 0000-0003-2059-5841
    affiliation: 1
  - name: Simon Stevenson
    orcid: 0000-0002-6100-537X
    affiliation: "4,5" 
  - name: Dorottya Szécsi
    orcid: 0000-0001-6473-7085
    affiliation: 6

    


affiliations:
 - name: Institute of Astronomy, KU Leuven, Celestijnenlaan 200D, B-3001, Leuven, Belgium
   index: 1
 - name: Department of Physics and Astronomy, University of North Carolina at Chapel Hill, 120 E. Cameron Avenue, Chapel Hill, NC 27599, USA
   index: 2
 - name: McWilliams Center for Cosmology, Department of Physics, Carnegie Mellon University, 5000 Forbes Avenue, Pittsburgh, PA 15213, USA
   index: 3
 - name: Centre for Astrophysics and Supercomputing, Swinburne University of Technology, Hawthorn, VIC 3122, Australia
   index: 4 
 - name: OzGrav-The ARC Centre of Excellence for Gravitational Wave Discovery, Hawthorn, VIC 3122, Australia
   index: 5
 - name: Institute of Astronomy, Faculty of Physics, Astronomy and Informatics, Nicolaus Copernicus University, Grudziadzka 5, 87-100, Torun, Poland
   index: 6


date: 16 Feb 2025
bibliography: paper.bib

---

# Summary

Stars, especially those greater than 9 solar masses (massive stars), play a pivotal role in shaping stellar populations. In the era of big-data astronomy, driven by the high-quality observational data produced by both ground-based and space-based electromagnetic telescopes, as well as gravitational wave and multi-messenger detectors, it is critical that the updates in stellar evolution can be easily incorporated into simulations that model their populations and interactions. Rapid stellar evolution codes offer a fast and computationally inexpensive way to capitalize on the population statistics provided by large surveys by integrating the latest theoretical models of stars into simulations of stellar populations.


# Statement of need

The current best way of computing stellar evolution is by solving equations of stellar structure in mass or radius coordinates through one-
Dimensional stellar codes (1D codes).  However, the inherent uncertainties in stellar evolution cause various 1D codes to adopt different physical inputs, leading to significant differences in the predictions for the evolution of stars and stellar populations. Moreover, computational requirements and robustness issues render 1D codes impractical for direct use in population synthesis codes. 

METISSE, written in Modern Fortran, interpolates between sets of pre-computed 1D stellar models to quickly calculate stellar parameters at each time step, thereby determining the stars' evolution for using stellar binary population synthesis codes. METISSE is fast and robust and the input stellar models can be easily swapped, allowing for systematic studies of stellar parameters on the stellar populations. It has already been used in a number of scientific publications, using stellar models from MESA and the Bonn code.  METISSE is similar to other interpolation-based rapid stellar evolution codes, such as SEVN [@Iorio et al. 2023] and Combine. However, METISSE is designed to serve as an alternative to the widely used SSE fitting formula [@Hurley:2000] and enables seamless integration with population synthesis codes that currently use SSE to calculate stellar parameters. 

METISSE can be used either as a standalone code (for single stellar populations) or in conjunction with population synthesis codes (for modelling binaries and star clusters). Additionally, multiple ongoing projects employ METISSE with Cbinary evolution codes with COSMIC, 
 With METISSE, one can combine diverse observations of stars and stellar populations into a single modelling framework capable of testing both binary and single stellar evolution physics. 


# Acknowledgements


# References