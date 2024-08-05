# Glossary

## Stellar phases

| Acronym | Stellar phase                         |
|-------|-------------------------------------------------|
| MS    | Main Sequence                                   |
| HG    | Hertzsprung Gap                                 |
| RGB   | Red Giant Branch                                |
| cHeB  | Core Helium Burning                             |
| AGB   | Asymptotic Giant Branch                         |
| TPAGB | Thermally Pulsing Asymptotic Giant Branch       |
| He-MS | Helium Main Sequence (stripped stars only)      |
| He-HG | Helium Hertzsprung Gap (stripped stars only)    |
| He-GB | Helium Giant Branch (stripped stars only)       |


## Output Files

METISSE can produce two types of output files:


### SSE-style file


SSE-style output files are text files ending with .dat and are controlled by the `write_output_to_file` function in SSE_input_controls. These files have a file structure similar to the popular [SSE](https://astronomy.swin.edu.au/~jhurley/bsedload.html) code. They contain the below listed stellar parameters up to the maximum age. Time and age at hydrogen ZAMS are assumed to be zero at the beginning of the code. 

::: {Important}
These files can only be created in standalone mode of METISSE.
:::

| Column Header | Description |
|-----------------|-----------------|
| time | Physical time [Myr] |
| age | Age of star [Myr] |
| mass | Current mass of the star [M<sub>$_\odot$</sub>] |
| core_mass | Mass of dominant core [M<sub>$_\odot$</sub>] |
| He_core | Mass of helium core [M<sub>$_\odot$</sub>] |
| CO_core | Mass of carbon-oxygen core [M<sub>$_\odot$</sub>] |
| log_L | Log of surface luminosity [L<sub>$_\odot$</sub>] |
| log_Teff | Log of effective temperature [K] |
| log_radius | Log of radius [R<sub>$_\odot$</sub>] |
| phase | SSE stellar type/phase |

### MIST-style file


MIST-style files, which end with .eep, are produced for debugging purposes. The file structure follows the same format used in the files from the [MIST](https://waps.cfa.harvard.edu/MIST/model_grids.html) project. METISSE can write a mass-interpolated track to an output file with the same columns as input files, including a stellar phase column. This output file only contains data from ZAMS to the end of nuclear-burning phases and does not include information about the remnant phase. It is controlled by the `write_eep_file` function in METISSE_input_controls.



## Metallicity File

This file ends with `_metallicity.in`and contains information about the input tracks for a specific metallicity, including the path to the folder containing the EEP tracks, the metallicity value, and other relevant information and metadata. For more details see [](using_custom_input_tracks.md#metallicity-controls).


## Format files


The format file contains information on how to read files, such as header location, important EEPs, and certain column names. The path to a format file is provided through the metallicity file. For more details see [](using_custom_input_tracks.md#format-controls).


## Fortran namelists

Namelists offer a convenient method for inputting data into Fortran programs. Each entry within a namelist follows this structure:

```
name = value ! comment
```

Values are provided using typical Fortran syntax. `SSE_input_controls` , `METISSE_input_controls`, `metallicity_controls` and `format_controls` are Fortran namelists, so comments (anything after the exclamation mark symbol!) and whitespaces can be used freely. Characters are also case-insensitive. It is important to ensure that there is a blank line at the end of the file following the `/` symbol.

:::{Note}

Never modify any file directly inside the `defaults` folder.

:::

## Stellar tracks and stellar models

A stellar model represents the physical parameters of a star at a specific point in time. By calculating a series of stellar models from the star's formation to the end of its nuclear-burning phase, we can determine the evolutionary history of the star, also known as its stellar track. In this documentation, we use "stellar track" to refer to a time-sequence of stellar models. In contrast "set of stellar models" or "set of stellar tracks" refers to the evolutionary tracks of stars with different initial masses but the same metallicity.



