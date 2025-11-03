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
    affiliation: 1 
  - name: Katelyn Breivik
    orcid: 0000-0001-5228-6598
    affiliation: 3 
  - name: Jarrod Hurley
    orcid: 0000-0003-2694-0415
    affiliation: "4,5" # (Multiple affiliations must be quoted)
  - name: Carl Rodriguez
    orcid: 0000-0003-4175-8881
    affiliation: 2 
  - name: Simon Stevenson
    orcid: 0000-0002-6100-537X
    affiliation: "4,5" 
  - name: Alex Kemp
    orcid: 0000-0003-2059-5841
    affiliation: 1
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


date: 3 Nov 2025
bibliography: paper.bib

---

# Summary

METISSE is an open-source stellar evolution tool specifically designed to be integrated with binary evolution and population synthesis codes. Unlike traditional rapid stellar evolution schemes based on the fitting formulae (e.g. SSE; [@Hurley:2000]), METISSE interpolates between pre-computed stellar models to quickly derive stellar parameters for population synthesis codes. This approach makes it straightforward to incorporate the latest stellar models and explore the effects of different input physics. METISSE can be easily integrated with any population synthesis code that currently uses popular Fortran code SSE to calculate stellar parameters. While METISSE can also be used in standalone mode to evolve single stellar populations, its primary purpose is to act as an interpolation engine within existing population synthesis codes for modelling binary stars and star clusters. Written in Modern Fortran, METISSE is both fast and robust, making it well-suited for large-scale population studies that require realistic stellar evolution input. The full documentation for METISSE is available at [metisse.readthedocs.io](https://metisse.readthedocs.io).

# Statement of need

Stars, especially those with masses greater than eight solar masses (massive stars), play a pivotal role in shaping stellar populations. The best way of computing stellar evolution involves solving equations of stellar structure and evolution through detailed stellar evolution codes such as MESA [@Jermyn:2023]. However, the inherent uncertainties in stellar evolution cause stellar codes to adopt different physical inputs, leading to significant differences in the predictions for the evolution of stars and stellar populations [@Agrawal:2022a]. Moreover, computational requirements and robustness issues render these codes impractical for direct use in large population synthesis simulations. 

Rapid stellar-evolution codes such as SSE [@Hurley:2000]that rely on fitting formulae manually calibrated to specific stellar models, have long provided a fast and efficient way to calculate stellar population properties. However, because the formulas must be recalculated manually for each new set of stellar models, they cannot be easily updated to incorporate advances in stellar evolution, limiting their flexibility and applicability to modern studies.

In recent years, interpolation-based rapid stellar-evolution codes such as TRILEGAL [@Girardi:2005], ComBinE [@Kruckow:2018], SEVN [@Iorio:2022], and MINT [@Rees:2025], along with machine-learning–based frameworks like POSYDON [@Fragos:2022], have become increasingly common as alternatives to SSE. However, these tools can not be easily incorporated into existing frameworks: binary population codes such as BSE [@Hurley:2002] and COSMIC [@Breivik:2020], or star cluster codes such as CMC [Rodriguez:2021] and NBODY6 [Aarseth:2003], that continue to rely on standard SSE routines.

METISSE provides performance and capabilities comparable to other interpolation-based codes but is implemented in Modern Fortran with full Fortran 77 compatibility. It includes SSE-style subroutines for straightforward integration into existing binary and population-synthesis codes. METISSE comes with an example set of [MESA stellar models](https://zenodo.org/records/17513335) and can also use other published grids, such as MIST [@Choi:2016] and BoOST [@Szecsi:2020]. Its interpolation framework allows stellar model grids to be easily swapped and explore how stellar parameters affect population outcomes. METISSE thus offers a flexible and efficient framework for stellar-evolution calculations in population synthesis, offering higher physical accuracy than analytic prescriptions without the heavy computational cost of full stellar-structure calculations. 

METISSE has already been employed in several scientific publications. For instance, it has also been used to demonstrate the impact of core overshooting — one of the major uncertainties in stellar evolution — on the evolutionary outcomes of binary systems [@Agrawal:2023]. Additionally, it has been used with stellar models from MESA as well as models from the Bonn Code (via the BoOST project [@Szecsi:2020]) to conduct a systematic study study of how different physical parameters affect the evolutionary properties of massive single stars [@Agrawal:2020]. 

Multiple ongoing projects use METISSE alongside the binary population synthesis code COSMIC[@Breivik:2020] to investigate the population properties of black hole-X-ray binaries, LISA white dwarf binaries, and GAIA black hole-star systems. 
In the era of big-data astronomy, driven by high-quality observational data from both ground-based and space-based telescopes, as well as gravitational wave and multi-messenger detectors, METISSE facilitates the seamless incorporation of updates in stellar evolution into simulations that model stellar populations and their interactions.

# Acknowledgements
We thank Duncan P. Maclean, Christopher Crow, Runqiu Ye and Steven Rieder for their help with testing METISSE. PA acknowledges support from the European Research Council (ERC) under the Horizon Europe programme (Synergy Grant agreement 101071505: 4D-STAR). While partially funded by the European Union, views and opinions expressed are however those of the author only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them. PA, JH and SS acknowledge support from the Australian Research Council Centre of Excellence for Gravitational Wave Discovery (OzGrav), through project number CE170100004. This research was supported by NSF grant AST-2310362. CR acknowledges support from the Alfred P. Sloan Foundation and the David and Lucile Packard Foundation. AK has received funding from the KU Leuven Research Council (grant C16/18/005: PARADISE). SS is supported by the ARC Discovery Early Career Research Award DE220100241. DSz acknowledges support from the National Science Center (NCN), Poland under grant number OPUS 2021/41/B/ST9/00757.


# References