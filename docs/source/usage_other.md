# Using METISSE with other codes


Presently, METISSE can be used with the following codes: 

1. BSE (Hurley et al. 2002)

To run METISSE with BSE, download the METISSE enabled version of BSE from [here](https://github.com/poojanagrawal/BSE-METISSE) and set `use_SSE = .false.` in the `bse.f` or `popbin.f`. 

The values of mass, metallicity and other binary parameters are read through one of the BSE input files. 
METISSE related inputs are read through the `METISSE_input_controls` inlist in the *evolve_metisse.in* file. 

METISSE looks for files(s) ending with *metallicity.in* in the paths provided through 'tracks_dir' and 'tracks_dir_he' where 'tracks_dir' is for providing the location of the metallicity file(s) for hydrogen star and  'tracks_dir_he' is for providing the location of the metallicity file(s) for the naked helium/stripped stars.

It checks for the following condition to find a match in metallicity 
`(abs(Z_input-Z_required)/MIN(Z_input,Z_required)) > Z_accuracy_limit`, 
 where `Z_input` is the metallicity value of the tracks (contained in the metallicity file) and `Z_required` is the value we want.

For a grid of stellar tracks, containing sets of stellar tracks with different metallicities, the above paths can contain a list of metallicity files, and METISSE will select the one that is closest to the input metallicity according to the aforementioned condition. 

If helium EEP tracks are not provided, METISSE will revert to using SSE formulae for the evolution of the naked helium stars.


2. COSMIC (Breivik et al. 2020) 

Running METISSE with COSMIC requires additional code from [COSMIC](https://github.com/COSMIC-PopSynth/COSMIC).


To use METISSE with COSMIC, set `stellar_engine` to `metisse` in SSEDict of COSMIC's input and supply paths for folders containing `metallicity files` for hydrogen and helium EEP tracks. 


`METISSE_input_controls` is not read when using COSMIC. Therefore `path_to_tracks` and `path_to_he_tracks` should also be provided in the SSEDict (They are used in place of the `tracks_dir`  and `tracks_dir_he` options of `METISSE_input_controls`). 


``` python
SSEDict = {'stellar_engine': 'metisse', 'path_to_tracks': path_to_tracks, 'path_to_he_tracks':path_to_he_tracks }
```

Here too, if `path_to_he_tracks` is empty, METISSE will revert to using SSE formulae for the evolution of the naked helium stars.

<!--## Example-->
