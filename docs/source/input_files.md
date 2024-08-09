# Input files

Inputs to METISSE are provided using the following namelists:


## SSE input controls

`SSE_input_controls` is holds input parameters describing the initial conditions of the star/stellar populations. For the most up-to-date variable names and their default values refer to [main_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/main_defaults.inc).

:::{Note}
`main.input` is **only** read in METISSE's standalone mode. 
When METISSE is used with other codes, the input parameters from the overlying code are used.
:::



### EVOLUTION CONTROLS

| Variable name        | Description  |
|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| read_mass_from_file  | If `.true.`, use initial masses listed in the `input_mass_file`.<br> Default is `.false.`.                                                                                                                                                |
| input_mass_file      | Path to the input mass file if `read_mass_from_file = .true.`. List one mass value per line of the file.                                                                                       |
| number_of_tracks     | Number of stars to be evolved. If `read_mass_from_file = .true.`, <br>it denotes the number of mass values to be read from the `input_mass_file`.                                                                                                                            |
| max_mass             | Upper limit of the unform distribution if `read_mass_from_file = .false.` and `number_of_tracks > 1`                                                                                                                                                               |
| min_mass             | Lower limit of the unform distribution if `read_mass_from_file = .false.`                                                                                              |
| initial_Z            | Initial metallicity                                                                                                                                        |
| max_age              | Maximum age in Myrs    |


:::{Note}

1. To evolve a single star set `number_of_tracks = 1` and `min_mass`, to the desired initial mass. In this case, `max_mass` is ignored.

2. To compute evolution of a population of stars with uniformly distributed initial masses, use `min_mass` to set the lower limit of the distribution, `max_mass` for the upper limit, and `number_of_tracks` for the size of the population. 
:::

```
! EVOLUTION CONTROLS
read_mass_from_file = .false.     
input_mass_file = ''
number_of_tracks = 0
max_mass = -1.0
min_mass = -1.0
initial_Z = -1.0
max_age = -1.0   
```

    
### REMNANT CONTROLS



| Variable name                 | Description |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| WD_mass_scheme                | White Dwarf (WD) luminosity calculation method:<br>(1) "Mestel" - [Shapiro S. L., Teukolsky S. A., 1983](https://ui.adsabs.harvard.edu/abs/1983bhwd.book.....S/abstract)<br>(2) "Modified_mestel" - [Hurley J. R., Shara M. M., 2003](https://iopscience.iop.org/article/10.1086/374637) <br> Default is "Modified_mestel"                                                                                              |
| use_initial_final_mass_relation | If `.true.` use the initial final mass relation for white dwarfs <br> from [Han, Z, Podsiadlowski, P., Eggleton, P. P. 1995](https://ui.adsabs.harvard.edu/abs/1995MNRAS.272..800H/abstract). <br> Default is `.false.`                                                                                               |
| BHNS_mass_scheme              | Neutron Star/Black Hole (NS/BH) type and mass calculation method:<br>(1) "original_SSE" - [Hurley et al. 2000](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract)<br>(2) "Belczynski2002" - [Belczynski et al. 2002](https://iopscience.iop.org/article/10.1086/340304)<br>(3) "Belczynski2008" - [Belczynski et al. 2008](https://iopscience.iop.org/article/10.1086/521026)<br>(4) "Eldridge_Tout2004" - [Eldridge J. J., Tout C. A., 2004](https://ui.adsabs.harvard.edu/abs/2004MNRAS.353...87E/abstract) <br> Default is "Belczynski2008"                                                                                            |
| max_NS_mass | Maximum neutron star mass. Recommended 1.8 for <br> BHNS_mass_scheme = "original_SSE", 3.0 otherwise. Default is 3.0                                                                                             |
| allow_electron_capture | Allow electron capture supernovae if `.true.`. <br> Default is `.true.`|        

```
! REMNANT CONTROLS 
WD_mass_scheme = 'Modified_mestel'
use_initial_final_mass_relation = .false.       
BHNS_mass_scheme = 'Belczynski2008'
max_NS_mass = 3.d0
allow_electron_capture = .true.  
```
### TIMESTEP CONTROLS

Like SSE, METISSE  uses `pts` variables to determine the timestep as the decimal fractions of the time spent in a phase. 

| Variable name                 | Description |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pts_1                         | Determine timestep for 95% of MS, HeMS. Default is 0.05
| pts_2                         | Determine timestep for last 5% of MS, cHeBurn, HeHG, HeGB. Default is 0.01
| pts_3                         | Determine timestep for HG, RGB, EAGB, TPAGB. Default is 0.02

```
! TIMESTEP CONTROLS
pts_1 = 0.05
pts_2 = 0.01
pts_3 = 0.02
```

### OUTPUT CONTROLS
| Variable name                 | Description |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| write_output_to_file           | Generate a [SSE-style](acronyms_definitions.md#sse-style-file) output file (file ending in .dat) at the END of the evolution. <br> Default is `.false.`.|

```
! OUTPUT CONTROLS
write_output_to_file = .false.
```

By default, METISSE only uses the essential columns and the additional columns specified in `format_file` for interpolation. However, when using METISSE in its standalone mode, there exists an option for interpolating in all the columns or certain extra columns of the input tracks. The interpolated quantities are printed in a [MIST-style](acronyms_definitions.md#mist-style-file) file if `write_eep_file = .true.` in `METISSE_controls`.


| Parameter            | Description     |
|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| read_all_columns     | If `.true.`, interpolate in all columns of input files. Default is `.false.`                                                                        |
| extra_columns        | Extra columns to be used for interpolation. Up to 100 column names can be listed <br> as strings, separated by commas. <br> *Only read if `read_all_columns` is `.false.`*                    |
| extra_columns_file   | Alternative way for listing extra column names. List one column name per line <br>in the `extra_columns_file` and provide path of the file here. <br> *Only read if `read_all_columns` is `.false.`*.                                                             |


```
read_all_columns = .false.
extra_columns = ''
extra_columns_file = ''

```


## METISSE input controls

`METISSE_input_controls` contains input parameters specific to METISSE. 

For the most up-to-date variable names and their default values refer to [metisse_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/metisse_defaults.inc). 


### TRACK CONTROLS


| Variable name   | Description  |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| METALLICITY_DIR     		| Absolute path to the folder containing the [metallicity files](acronyms_definitions.md#metallicity-file) for the hydrogen stars. <br> **METISSE will stop if `METALLICITY_DIR` is an empty string.**                                                                                                                                      |
| METALLICITY_DIR_HE  		| Absolute path to the folder containing the [metallicity files](acronyms_definitions.md#metallicity-file) for the naked helium stars. <br>If `METALLICITY_DIR_HE` is an empty string, SSE fitting formulae for naked helium stars <br> are used.                                                                                                                                      |
| Z_accuracy_limit		| METISSE checks for a match in metallicity based on the condition: <br> $ \frac{abs(Z\_input - Z\_required)}{min(Z\_input, Z\_required)} >$ `Z_accuracy_limit`, <br> where $Z\_input$ is the metallicity value of the tracks (from the metallicity file)<br>   and $Z\_required$ is the desired value. Default is 10<sup>-2</sup>                                                                                                       |
| mass_accuracy_limit  	| METISSE skips mass-interpolation and uses a neighboring mass track if <br> the absoulte difference between the initial mass of that track <br>and the input mass is less than the `mass_accuracy_limit`. Default is 10<sup>-4</sup>|


```
! TRACK CONTROLS
METALLICITY_DIR = ''
METALLICITY_DIR_HE = ''
Z_accuracy_limit = 1d-2
mass_accuracy_limit = 1d-4

```

    
### MISCELLANEOUS CONTROLS
| Variable name                 | Description  |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| verbose                       | If `.true.`, print important details when reading the input tracks on the terminal. <br>Default is `.false.`.| 
| write_eep_file                | Generate a [MIST-style](acronyms_definitions.md#mist-style-file) output file at every step of mass interpolation. <br> Useful when computing evolution of isolated single stars with <br>implicit mass loss and for debugging purposes. Default is `.false.`.|
| write_error_to_file           | If `.true.`, error messages from METISSE are written to a file named `fort.99` <br>otherwise they are printed on screen (does NOT apply to errors during <br>reading input files). Default is `.true.`.|
| construct_postagb_track            | Artificially construct the HRD (Hertzsprung-Russell diagram) track between <br> tip of the AGB to the WD cooling track (for low-mass stars). It is <br>important if input tracks do not contain this phase but can be used regardless. <br>Default is `.false.`.| 



```
! MISCELLANEOUS CONTROLS

verbose = .false.
write_eep_file = .false.		
write_error_to_file = .true.
construct_postagb_track = .false.

```



