# Using custom input tracks


METISSE can use any set of stellar tracks computed with different stellar evolution codes. The only requirement is that these tracks should converted to the equivalent evolutionary point (EEP) format before use in METISSE. In EEP-format, significant evolutionary points such as the zero-age main sequence (ZAMS) or terminal age main sequence (TAMS) occur at the same line number across each file. 
Stellar tracks can be easily converted to EEP format using code packages such as [ISO](https://github.com/aarondotter/iso). Important details about these tracks, such as the their metallicity value, file structure, names of certain major columns should also be provided through the inlists `metallicity_controls` and `format_controls`.


<!-- :::{Important}

While different sets of EEP tracks can share the same format file, they **must** have separate metallicity files.
::: 
 -->


## Metallicity controls 

Details pertaining to the metallicity value of the EEP tracks and location of the files are provided through the `metallicity_controls` inlist, found in the relevant `metallicity file`. 

:::{Note}

1. Different sets of EEP tracks **must** have separate metallicity files.

2. The name of a `metallicity file` can be arbitary but must end with `_metallicity.in` (for example, filename_metallicity.in).

::: 

For the most up-to-date variable names and their default values refer to [metallicity_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/metallicity_defaults.inc).


The following three parameters are essential for any given set of input tracks. 
**METISSE will stop 
and raise an error if either of these is not provided.**   


| Parameter            | Description     |
|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| INPUT_FILES_DIR      | Location of the folder containing the EEP tracks relative to the metallicity file.                                                                             |
| Z_files              | Metallicity value of the EEP tracks. This is cross-matched against the input metallicity value.                               |
| format_file          | Location of the format file relative to the metallicity file. |


```

INPUT_FILES_DIR = ''
Z_files = -1.0
format_file = ''
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

| Parameter            | Description     |
|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Z_H                  | Hydrogen abundance. Default is SSE formulae. <br> If < 0, calculated from Z as $0.76-3\times Z$.                                                                                |
| Z_He                 | Helium abundance. Default is SSE formulae. <br>If < 0, calculated from Z as $0.24+2\times Z$. |

```
Z_H = -1.0
Z_He = -1.0

```

 `metallicity_controls` also contains the option to provide user-defined values of mass cutoff parameters or zparameters (See [Agrawal et al. 2020](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..933A/abstract) for details). If negative, these values are calculated by the code. 

| Parameter            | Description     |
|----------------------|----------------------------------------------------------|
| Mhook                | Mass below which hook does not appear on the MS.			|
| Mhef                 | Mass above which helium ignition occurs non-degenerately.              |
| Mfgb                 | Mass above which helium ignition occurs on the HG.                     |
| Mup                  | Mass below which carbon ignition does not occur.                    |
| Mec                  | Mass above which carbon ignites in the centre. |                 
| Mextra               | Extra mass cutoff (if any).  |

```
Mhook = -1.0
Mhef = -1.0
Mfgb = -1.0
Mup = -1.0
Mec = -1.0
Mextra = -1.0   
```


## Format controls

Apart from the metallicity files, users are also required to specify the structure/format of the files conatining input EEP tracks through the `&format_controls` inlist. For the most up-to-date variable names and their default values refer to [format_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/format_defaults.inc). 

| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| read_eep_files          | Set to `.true.` only for [MIST-style](https://waps.cfa.harvard.edu/MIST/model_grids.html) files; set to `.false.` for other types of files. <br> If `.true.`, other options in this section do not need to be provided. The required details are read directly from the input files <br> .                                                                                                                                         |
| file_extension          | File extension of the input files (e.g., '.eep', '.dat'). METISSE will look for files ending with `file_extension` in the `INPUT_FILES_DIR`.                                                                                          |

| header_location         | Line number corresponding to the column names in the input files for EEP-tracks. <br> If the input files do not contain column names; specify `column_name_file` for such cases.                    																	 |
| extra_char              | Any extra character present at the beginning of the header line.                                                                                                         |
| column_name_file        | File containing column names (one per line) if `header_location` <= 0 <br> and column names cannot be determined from the input files.                                                                                             |
| total_cols              | Total number of columns in the input files.                            |

```

file_extension = ''
read_eep_files = .false.
header_location = -1
extra_char = ''
column_name_file = ''
total_cols = -1

```


    
### Essential columns


Use this section to specify essential column names. In the case of the luminosity and the radius columns, at least one of the absolute or log column names should be provided. 
Make sure that the units are correct.
**METISSE will stop and raise an error if these columns can not be located.**



| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| age_colname             | Column name for age in years.                                                                       |
| mass_colname            | Column name for star's total mass in solar units.                                                                                                                               |
| log_L_colname           | Column name for log10 of luminosity in solar units.                                                                                    |
| Lum_colname             | Column name for stellar luminosity in solar units. <br>*Used only if `log_L_colname` is not provided.*                                                                  |
| log_R_colname           | Column name for log10 of radius in solar units.                                                           |
| Radius_colname          | Column name for stellar radius in solar units. <br> *Used only if `log_R_colname` is not provided.*                                                                   |
| he_core_mass            | Column name for mass of He enriched/H depleted core in solar units.                                                                                                            |
| co_core_mass            | Column name for mass of C enriched/He depleted core in solar units.                               |



```
age_colname = ''
mass_colname = ''
log_L_colname = ''
Lum_colname = ''
log_R_colname = ''
Radius_colname = ''
he_core_mass = ''
co_core_mass = ''

```

### Additional columns 

METISSE will raise error but not stop if these column names are not provided. However, it will revert to using SSE formulae where needed. 



| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| he_core_radius          | Column name for radius of helium enriched/hydrogen depleted core <br>in solar units (cannot use log).                                                                                                                                                                      |
| co_core_radius          | Column name for radius of carbon enriched/helium depleted core <br>in solar units (cannot use log).                                                                         |
| mass_conv_envelope      | Column name for mass of the convective envelope in solar units.                                                                                               |
| radius_conv_envelope    | Column name for radius of the convective envelope in solar units.                                                                                              |
| log_T_colname           | Column name for log10 of surface temperature in K.                                                                    |
| Teff_colname            | Column name for surface temperature in K. <br> *Used only if `log_T_colname` is not provided.*                                                                      |
| log_Tc                  | Column name for central temperature in log units.                                                                                                            |
| he4_mass_frac           | Column name for helium-4 mass fraction at centre.                                                                                                            |
| c12_mass_frac           | Column name for carbon-12 mass fraction at centre.                                                                                                           |
| o16_mass_frac           | Column name for oxygen-16 mass fraction at centre.     |


```
he_core_radius = ''    
co_core_radius = ''
mass_conv_envelope = ''
radius_conv_envelope = ''
log_T_colname = ''
Teff_colname = ''
log_Tc = ''
he4_mass_frac = ''
c12_mass_frac = ''
o16_mass_frac = ''
```
### EEP details for hydrogen stars

From a set of input models, METISSE needs to know the locations of key EEPs in order to assign stellar phases to the interpolated tracks. These phases are identical to those used by the SSE code (See Table 1. of [Agrawal et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.4549A/abstract)) and are important for certain decision-making processes in the code, particularly for binary evolution.

| EEP    			      | Corresponding evolutionary point     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| PreMS_EEP               | Pre-Main Sequence.   			 |
| ZAMS_EEP                | Zero-Age Main Sequence.                                                                                                                  |
| IAMS_EEP                | Intermediate-Age Main Sequence.                                                                                                                 |
| TAMS_EEP                | End of the Main Sequence.                                                                                                                 |
| BGB_EEP                 | Base of the Giant Branch.                         |
| cHeIgnition_EEP         | Helium Ignition in the core.                                                                                                                  |
| cHeBurn_EEP             | Core Helium Burning.                                                                                                                   |
| TA_cHeB_EEP             | End of core Helium Burning.                                                                                                |
| TPAGB_EEP               | Thermally Pulsing Asymptotic Giant Branch.   |
| cCBurn_EEP              | End of core Carbon Burning.                                                                                                                    |
| post_AGB_EEP            | Post-Asymptotic Giant Branch.      |                                                                    |
| Extra_EEP1              | Additional EEP *(optional)*.                                                                                                                        |
| Extra_EEP2              | Additional EEP *(optional)*.                                                                                                                        |
| Extra_EEP3              | Additional EEP *(optional)*.                                                                                                                          |
| Initial_EEP             | Line number to start reading files from. <br> If < 0, ZAMS_EEP is used.                                                                |
| Final_EEP               | Line number to stop reading files at. <br> If < 0, maximum of the listed EEPs is used.                                                                                                                     


```
PreMS_EEP = -1
ZAMS_EEP = -1
IAMS_EEP = -1
TAMS_EEP = -1
BGB_EEP = -1
cHeIgnition_EEP = -1
cHeBurn_EEP = -1
TA_cHeB_EEP = -1
TPAGB_EEP = -1
cCBurn_EEP = -1
post_AGB_EEP = -1

Extra_EEP1 = -1
Extra_EEP2= -1
Extra_EEP3 = -1

Initial_EEP = -1
Final_EEP = -1

```


### EEP details for stripped/ naked helium stars

Locations of key EEPs for stars that have lost their hydrogen-rich envelopes, also known as naked helium or stripped stars. Similar to primary EEPs of hydrogen stars, these are used to assign stellar phases to the interpolated tracks for naked helium star phases.


| EEP     				  | Corresponding evolutionary point     |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| ZAMS_HE_EEP             | Zero Age Main Sequence of Helium stars.                                                                                               |
| TAMS_HE_EEP             | Terminal Age Main Sequence of Helium stars.                                                                                             |
| GB_HE_EEP               | Base of the Giant Branch of Helium stars.                                                                                                   |
| cCBurn_HE_EEP           | end of core Carbon Burning of Helium stars .                                                                                                 |
| post_AGB_HE_EEP         | Post-Asymptotic Giant Branch of Helium stars.                                                                                            |

```
ZAMS_HE_EEP = -1
TAMS_HE_EEP = -1
GB_HE_EEP = -1
cCBurn_HE_EEP = -1
post_AGB_HE_EEP = -1 
```


### How to deal with the incomplete tracks

:::
> *"In an ideal world, we wouldn't need this section, but unfortunately, the world is not ideal, and incomplete or incorrect data is more common in datasets than we'd like to believe."*
:::

For a given initial mass, METISSE calculates an evolutionary track by interpolating between the corresponding EEPs of neighbouring mass tracks. 
If any of the neighbouring tracks is incomplete, the interpolated track is also rendered incomplete. This section explains how METIISE identifies incomplete tracks and attempts to fill in the missing data. 

| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| low_mass_final_eep      | Final EEP for stars with M < Mec .                                                                    |
| high_mass_final_eep     | Final EEP for stars with M >= Mec .                                                                      |
| fix_track               | If `.true.`, METISSE relaxes criteria for finding neighboring tracks to fix <br> incomplete tracks.                                                                     |
| lookup_index            | Determines the mass range for searching neighboring tracks when fixing <br>incomplete tracks. The range is M-(M$\times$`lookup_index`) and M+(M$\times$`lookup_index`), <br>where M is the initial mass of the star. |


```    

low_mass_final_eep = -1
high_mass_final_eep = -1
fix_track = .true.
lookup_index = 1.0

```