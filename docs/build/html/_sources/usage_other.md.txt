# Using METISSE to evolve binaries

**METISSE**, by itself, is not capable of evolving stars in binary systems.  
However, it has been specifically designed to be integrated with binary evolution codes that perform such calculations, allowing them to make use of modern stellar models in binary evolution studies.

In particular, METISSE can be seamlessly integrated with codes that currently use the **SSE fitting formulae** [(Hurley et al. 2000)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) for single-star evolution computations.

Currently, METISSE is available as a module within the following codes:

```{toctree}
bse.md
cosmic.md
