# Using custom input tracks


METISSE can use any set of stellar tracks computed with different stellar evolution codes. The only requirement is that these tracks should converted to the equivalent evolutionary point (EEP) format before use in METISSE. In EEP format, significant evolutionary points such as the zero-age main sequence (ZAMS) or terminal-age main sequence (TAMS) occur at the same line number across each file. 
Stellar tracks can be easily converted to EEP format using code packages such as [ISO](https://github.com/aarondotter/iso). Important details about these tracks, such as their metallicity value, file structure, and names of certain major columns should also be provided through the inlists `metallicity_controls` and `format_controls`.


## Metallicity controls 

Details pertaining to the metallicity value of the EEP tracks and location of the files are provided through the `metallicity_controls` inlist, found in the relevant *metallicity file*. 

:::{Note}

1. Different sets of tracks **must have different metallicity files**, even if the format file referenced by each metallicity file is the same.

2. The name of a `metallicity file` can be arbitrary but must end with `_metallicity.in` (for example, filename_metallicity.in).

::: 

For the most up-to-date variable names and their default values refer to [metallicity_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/metallicity_defaults.inc).


The following three parameters are essential for any given set of input tracks. 
**METISSE will stop and raise an error if either of these is not provided.**   


| Parameter            | Description     |
|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| eep_tracks_dir     | Location of the folder containing the EEP tracks relative to the metallicity file.                                                                             |
| Z_files              | Metallicity value of the EEP tracks. This is cross-matched against the input metallicity value.                               |
| format_file          | Location of the format file relative to the metallicity file. |


```

eep_tracks_dir = ''
Z_files = -1.0
format_file = ''
```

| Parameter            | Description     |
|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Y_files              | Helium abundance. Default is SSE formulae. <br>If < 0, calculated from Z as $0.24+2\times Z$. |

```

Y_files = -1.0

```


METISSE identifies specific critical masses, referred to as Z-parameters or mass cutoffs. These critical masses mark the point where certain physical properties begin to manifest in stellar tracks (refer to appendix A of [Agrawal et al. 2020](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..933A/abstract) for more details).
Z-parameters are fixed for a given metallicity and can be determined through input stellar tracks using the essential columns and the additional columns supplied in the format file. 
However, if certain columns are missing from the input tracks or if the Z-parameters are not being located correctly, the user can input predetermined values for the Z-parameters.

| Parameter            | Corresponding initial mass of the star    |
|----------------------|----------------------------------------------------------|
| Mhook                | Above which hook feature starts to appear on the MS.			|
| Mhef                 | Above which helium ignition occurs non-degenerately in the core.              |
| Mfgb                 | Above which helium ignition occurs on the HG.                     |
| Mup                  | Above which off-centre C/O ignition can occur non-degenerately in the core.                    |
| Mec                  | Above which a star avoids electron captures on neon and proceeds to form an iron core. |                 
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

Apart from the metallicity files, users should also specify the format of the files containing input EEP tracks through the `&format_controls` inlist contained in what we call a *format file*. 
While different sets of EEP tracks can share the same format file, users should provide at least two format files: one for usual hydrogen stars and the other for helium star tracks.

For the most up-to-date variable names and their default values refer to [format_defaults](https://github.com/TeamMETISSE/METISSE/blob/develop/src/defaults/format_defaults.inc). 




### File details

| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| read_eep_files          | Set to `.true.` only for [MIST-style](https://waps.cfa.harvard.edu/MIST/model_grids.html) files; set to `.false.` for other types of files. <br> If `.true.`, other options in this section (until `total_cols`) do not need to be provided. <br> The required details are read directly from the input files.                                                                                                                                         |
| file_extension          | File extension of the input files (e.g., '.eep', '.dat'). <br> METISSE will look for files ending with `file_extension` in the `eep_tracks_dir`.                                                                                          |
| header_location         | Line number corresponding to the column names in the input files for EEP tracks. <br> If the input files do not contain column names; specify `column_name_file` for such cases.                    																	 |
| extra_char              | Any extra character present at the beginning of the header line.                                                                                                         |
| column_name_file        | File containing column names (one per line) if `header_location` <= 0 <br> and column names cannot be determined from the input files.                                                                                             |
| total_cols              | Total number of columns in the input files.                            |

```

read_eep_files = .false.
file_extension = ''
header_location = -1
extra_char = ''
column_name_file = ''
total_cols = -1

```


    
### Essential columns

As columns can be named differently by different stellar codes, METISSE needs to know certain important column names. In the case of the luminosity and the radius columns, at least one of the absolute or log column names should be provided. 
Make sure that the units are correct.
**METISSE will stop and raise an error if these columns can not be located.**



| Parameter               | Corresponding column name  |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| age_colname             | Age of the star in years.                                                                       |
| mass_colname            | Total mass of the star in solar units.                                                                                                                               |
| log_L_colname           | Log10 of stellar luminosity in solar units.                                                                                    |
| Lum_colname             | Stellar luminosity in solar units. <br>*Used only if `log_L_colname` is not provided.*                                                                  |
| log_R_colname           | Log10 of stellar radius in solar units.                                                           |
| Radius_colname          | Stellar radius in solar units. <br> *Used only if `log_R_colname` is not provided.*                                                                   |
| he_core_mass            | Mass of the helium enriched/hydrogen depleted core in solar units.                                                                                                            |
| co_core_mass            | Mass of the carbon enriched/helium depleted core in solar units.                               |


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

The following columns are useful in binary evolution calculations. 
METISSE will raise an error but not stop if these column names are not provided. However, it will revert to using SSE formulae where needed. 



| Parameter               | Corresponding column name    |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| he_core_radius          | Radius of the helium enriched/hydrogen depleted core <br>in solar units (cannot use log).                                                                                                                                                                      |
| co_core_radius          | Radius of the carbon enriched/helium depleted core <br>in solar units (cannot use log).                                                                         |
| mass_conv_envelope      | Mass of the convective envelope in solar units.                                                                                               |
| radius_conv_envelope    | Radius of the convective envelope in solar units.                                                                                              |
|binding_energy_colname	  | Binding energy of the envelope in ergs (cannot use log)		|


```
he_core_radius = ''    
co_core_radius = ''
mass_conv_envelope = ''
radius_conv_envelope = ''
binding_energy_colname = ''

```

The following columns are used in determining various Z-parameters.
If these are not provided, then user-defined values for Z-parameters (through the `metallicity_controls`) are used. If those are not provided either, then Z-parameters are determined using SSE formulae.


| Parameter               | Corresponding column name    |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| log_T_colname           | Log10 of surface temperature in K.                                                                    |
| Teff_colname            | Surface temperature in K. <br> *Used only if `log_T_colname` is not provided.*                                                                      |
| log_Tc                  | Central temperature in log units.                                                                                                            |
| he4_mass_frac           | Helium-4 mass fraction at centre.                                                                                                            |
| c12_mass_frac           | Carbon-12 mass fraction at centre.                                                                                                           |
| o16_mass_frac           | Oxygen-16 mass fraction at centre.     |


```

log_T_colname = ''
Teff_colname = ''
log_Tc = ''
he4_mass_frac = ''
c12_mass_frac = ''
o16_mass_frac = ''
```
### EEP details 
the 
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
| Initial_EEP             | EEP to start reading files from. <br> If < 0, ZAMS_EEP is used.                                                                |
| Final_EEP               | EEP to stop reading files at. <br> If < 0, the maximum of the listed EEPs is used.    |                                                                                                                 


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


When working with naked helium star tracks, all the above-mentioned format file options can be used to specify the file and the column details. However, only the following EEPs are used to assign naked helium star phases to the interpolated tracks: 


| EEP     				  | Corresponding evolutionary point for naked helium stars    |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| cHeBurn_EEP             | Zero Age Helium Main Sequence                                                                                               |
| TA_cHeB_EEP             | Terminal Age Helium Main Sequence |                                                                                           |
| BGB_EEP                 | Base of the Giant Branch                                                                                                   |
| cCBurn_EEP           	  | End of core Carbon Burning                                       
| TPAGB_EEP               | End of AGB |                                                           |
| post_AGB_EEP         	  | Post-AGB (end of AGB to WD cooling track) |
| Initial_EEP             | EEP to start reading files from. <br> If < 0, cHeBurn_EEP is used.                                                                |
| Final_EEP               | EEP to stop reading files at. <br> If < 0, the maximum of the listed EEPs is used.    |                                                                                                                 


Once read, these EEPs get stored with the appropriate (helium files specific) variable name inside METISSE. Other EEPs (such as ZAMS_EEP and TAMS_EEP) are ignored as they don't have an equivalent in helium star evolution.


### How to deal with incomplete tracks

:::
> *"In an ideal world, we wouldn't need this section, but unfortunately, the world is not ideal, and incomplete or incorrect data is more common in datasets than we'd like to believe."*
:::

For a given initial mass, METISSE calculates an evolutionary track by interpolating between the corresponding EEPs of neighbouring mass tracks. 
If any of the neighbouring tracks is incomplete, the interpolated track is also rendered incomplete. This section explains how METISSE identifies incomplete tracks and attempts to fill in the missing data. 

| Parameter               | Description     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| low_mass_final_eep      | Length of track for stars with initial mass < Mec is compared <br> against `low_mass_final_eep` to determine completeness.                                                                    |
| high_mass_final_eep     | Length of track for stars with initial mass >= Mec is compared <br> against `high_mass_final_eep` to determine completeness.                                                                      |
| fix_track               | If `.true.`, METISSE relaxes criteria for finding neighboring tracks to fix <br> incomplete tracks.                                                                     |
| lookup_index            | Determines the mass range for searching neighboring tracks when fixing <br>incomplete tracks. The range is M-(M$\times$`lookup_index`) and M+(M$\times$`lookup_index`), <br>where M is the initial mass of the star. |


```    

low_mass_final_eep = -1
high_mass_final_eep = -1
fix_track = .true.
lookup_index = 1.0

```

When evaluating the completeness of input tracks for naked helium stars, only the parameters `low_mass_final_eep` and `high_mass_final_eep` are considered and saved internally under helium-specific names. 
The options `fix_track` and `lookup_index` are not used in the format files of helium stars. Instead, the values provided in the format file for hydrogen stars are applied to both hydrogen and helium star tracks.

