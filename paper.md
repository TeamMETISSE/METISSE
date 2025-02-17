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

METISSE is an open-source code that interpolates between pre-computed stellar models to quickly derive stellar parameters for population synthesis codes. Written in Modern Fortran, METISSE is both fast and robust. It is comparable to existing rapid stellar evolution codes, such as SEVN [@Iorio:2023] and Combine [@Kruckow:2018] but with the added advantage that it can be seamlessly integrated with any population synthesis code that currently uses popular SSE fitting formulae[@Hurley:2000] to calculate stellar parameters. METISSE can function as a standalone code for simulating single stellar populations as well as in conjunction with population synthesis codes for modelling binary stars and star clusters. 

# Statement of need

Stars, especially those with masses greater than eight solar masses (massive stars), play a pivotal role in shaping stellar populations. The best way of computing stellar evolution involves solving equations of stellar structure and evolution through detailed stellar evolution codes such as MESA. However, the inherent uncertainties in stellar evolution cause stellar codes to adopt different physical inputs, leading to significant differences in the predictions for the evolution of stars and stellar populations[@Agrawal:2022a]. Moreover, computational requirements and robustness issues render these codes impractical for direct use in large population synthesis simulations. 

As a result, rapid stellar evolution codes that rely on fitting formulas, manually calibrated to specific sets of stellar models to determine the evolution of individual star systems, have been a popular alternative in the past. These methods provide a computationally efficient, fast, and reliable way to calculate stellar population properties. Unfortunately, fitting formulas must be recalculated manually each time for different stellar model sets, limiting our ability to conduct systematic studies of stellar parameters. The use of interpolation in METISSE allows input stellar models to be easily swapped and enables systematic studies of stellar parameters on the stellar populations. 

METISSE has already been employed in several scientific publications. For instance, it has been used to demonstrate the impact of core overshooting — one of the major uncertainties in stellar evolution — on the evolutionary outcomes of massive binary systems[@Agrawal:2023]. Multiple ongoing projects use METISSE alongside the binary population synthesis code COSMIC[@Breivik:2020] to investigate the population properties of black hole-X-ray binaries, LISA white dwarf binaries, and GAIA black hole-star systems. In the era of big-data astronomy, driven by high-quality observational data from both ground-based and space-based telescopes, as well as gravitational wave and multi-messenger detectors, METISSE facilitates the seamless incorporation of updates in stellar evolution into simulations that model stellar populations and their interactions.

# Acknowledgements


# References