# FAQ


## How to cite METISSE?
Please cite the following papers if you are using METISSE in your work. 


<!-- - TBA

More details about METISSE including code capabilities are described in the following papers: -->


- [Modelling stellar evolution in mass-transferring binaries and gravitational-wave progenitors with METISSE](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..933A/abstract)
- [The fates of massive stars: exploring uncertainties in stellar evolution with METISSE](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.4549A/abstract)

## Which interpolation method is used by METISSE? 

METISSE uses monotonic interpolation with a piece-wise cubic from [Steffen 1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S/abstract) to estimate stellar properties between pre-computed tracks. The method ensures smooth and accurate results but requires two points below and two points above the target value; consequently, METISSE switches to linear interpolation near the boundaries where sufficient points are not available.

<!-- The interpolation method is is fixed; changing it requires modifying the source code. -->

## Can I use METISSE to restart the evolution of binaries?
We are working hard to add this capability, but currently METISSE cannot restart the evolution of binaries with any of the supported overlying binary codes.


## Can METISSE interpolate in metallicity?
METISSE delibertely does not support interpolation in metallicity. However, you can choose the nearest metallicity value by adjusting `Z_accuracy_limit`. Please refer to [](usage_standalone.md#if-the-metallicity-value-is-not-present-in-the-input-grid) for more details.

## Can I add METISSE to my code?
Yes, if it uses Fortran based fitting formulae from SSE [(Hurley et al. 2000)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) to compute stellar parameters. Your code must contain Fortran77 subroutines - namely zcnsts.f, star.f, hrdiag.f, deltat.f, mlwind.f that come with the SSE code. For details refer to [](structure.md#adding-metisse-to-your-code).

## How can I contribute to METISSE?
We are always working on improving METISSE. 
If you would like to contribute by adding a new feature or creating with a bug fix, please 
1. [Fork](https://github.com/TeamMETISSE/METISSE/fork) the METISSE repository from GitHub. 
2. Make your changes in a feature branch.
3. Create a [pull request](https://github.com/TeamMETISSE/METISSE/pulls) with your proposed changes.

## I seek support or have to report problems with METISSE
For bug reports or issues, create an [issue](https://github.com/TeamMETISSE/METISSE/issues/new) on GitHub.
For general questions or discussion, use GitHub [discussions](https://github.com/TeamMETISSE/METISSE/discussions).