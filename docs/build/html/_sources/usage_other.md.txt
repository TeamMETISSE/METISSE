# Using METISSE with other codes

**Prequisite:** METISSE requires gcc/6.4.0 or above


Currently, METISSE can be used with the following codes: 

1. BSE (Hurley et al. 2002)

To run METISSE with BSE, simply enable `use_SSE = .false.` in the `bse.f` or `popbin.f`. 

The `metallicity_file_list` option in the `METISSE_input_controls` is used to specify the location of the metallicity file when METISSE is used in standalone mode or with BSE. Users can provide a list of metallicity files, and METISSE will select the one that is closest in input metallicity.
Similarly, EEP tracks for naked helium stars or stripped stars can be provided using `metallicity_file_list`. If helium EEP tracks are not provided, METISSE will revert to using SSE formulae for the evolution of the naked helium stars.

Other details are read through `METISSE_input_controls` inlist in the *evolve_metisse.in* file.


2. COSMIC (Breivik et al. 2020) 

Running METISSE with COSMIC requires additional code from [COSMIC](https://github.com/COSMIC-PopSynth/COSMIC).


To use METISSE with COSMIC, set `stellar_engine` to `metisse` in SSEDict of COSMIC's input and supply paths for folders containing `metallicity files` for hydrogen and helium EEP tracks. 


`METISSE_input_controls`  is not read when using COSMIC, therefore `metallicity_file_list` cannot be not directly provided. Instead, the `path_to_tracks` and `path_to_he_tracks` are specified SSEDict. METISSE then searches for all files ending with 'metallicity.in' in that location and creates the `metallicity_file_list`.


``` python
SSEDict = {'stellar_engine': 'metisse', 'path_to_tracks': path_to_tracks, 'path_to_he_tracks':path_to_he_tracks }
```

Here too, if `path_to_he_tracks` is empty, METISSE will revert to using SSE formulae for the evolution of the naked helium stars.

## Example
