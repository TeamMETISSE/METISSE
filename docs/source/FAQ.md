# FAQ

## How can I contribute to METISSE?
We are always working on improving METISSE. 
If you would like to contribute by adding a new feature or creating with a bug fix, please 
1. [fork](https://github.com/TeamMETISSE/METISSE/fork) the METISSE repository from GitHub. 
2. Make your changes in a feature branch.
3. Create a [pull request](https://github.com/TeamMETISSE/METISSE/pulls) with your proposed changes.

## I seek support or I have to report problems with METISSE
For bug reports or issues, create an [issue](https://github.com/TeamMETISSE/METISSE/issues/new) on GitHub.
For general questions or discussion, use GitHub [discussions](https://github.com/TeamMETISSE/METISSE/discussions).


## How to cite METISSE
Please cite the following papers if you are using METISSE in your work. 


<!-- - TBA

More details about METISSE including code capabilities are described in the following papers: -->


- [Modelling stellar evolution in mass-transferring binaries and gravitational-wave progenitors with METISSE](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..933A/abstract)
- [The fates of massive stars: exploring uncertainties in stellar evolution with METISSE](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.4549A/abstract)


## Can I use METISSE to restart the evolution of binaries?
We are working hard to add this capability, but currently METISSE cannot restart the evolution of binaries with any of the supported overlying binary codes.


## Can METISSE interpolate in metallicity?
METISSE does not support interpolation in metallicity. However, you can choose the nearest metallicity value by adjusting `Z_accuracy_limit`. Please refer to [](usage_standalone.md#if-the-metallicity-value-is-not-present-in-the-input-grid) for more details.


## How can I add METISSE to my code?
Adding METISSE to a code that currently uses SSE is very easy. Please get in touch through GitHub [discussions](https://github.com/TeamMETISSE/METISSE/discussions).


