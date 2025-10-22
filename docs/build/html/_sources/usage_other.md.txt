# Using METISSE to evolve binaries

METISSE, by itself, is not capable of evolving stars in binary systems. However, it has been specifically designed to be integrated with binary evolution codes that can do so, enabling them to make use of modern stellar data for binary evolution computations.

It can be seamlessly integrated with codes that currently use SSE fitting formulae [(Hurley et al. 2000)](https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract) for single star evolution computations.
It has already been added to the binary evolution codes BSE [(Hurley et al. 2002)](https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract) and COSMIC [(Breivik et al. 2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract).



## COSMIC 


1. To use METISSE with COSMIC one needs the GitHub repository for the latest version of COSMIC.

If you don't have COSMIC installed already you can fork it from the [GitHub repository](https://github.com/COSMIC-PopSynth/COSMIC/tree/develop) and clone it on your machine.


``` console

$ git clone https://github.com/COSMIC-PopSynth/COSMIC

```

Follow the COSMIC's documentation to [load the cosmic environment](https://cosmic-popsynth.github.io/COSMIC/pages/install.html)

Next, cd into the COSMIC's directory and compile the library: 


``` console

$ pip install .

```


Try evolving an example binary following the instructions on [this page](https://cosmic-popsynth.github.io/COSMIC/pages/examples.html). 


2. Next checkout the METISSE-Integrate branch of COSMIC. In the COSMIC's directory, do:


``` console

$ git checkout METISSE-integrate

```

Get the METISSE submodule.

``` console

$ git submodule update --init --recursive

```


Re-compile the library


``` console

$ pip install .

```

3. METISSE is now ready to use with COSMIC. 

The usual input namelist `METISSE_input_controls` is not read when using COSMIC. Instead a Python dictionary, `SSEDict` is used to set the `stellar_engine` and provide the location of the folders containing [metallicity files](glossary.md#metallicity-file) for hydrogen and helium stars through variables `path_to_tracks` and `path_to_he_tracks`.

``` ipython
SSEDict = {
          'stellar_engine': 'metisse', 
          'path_to_tracks': '/Users/poojan/Downloads/sample_tracks_solarZ/Hydrogen/', 
          'path_to_he_tracks': '/Users/poojan/Downloads/sample_tracks_solarZ/Helium/' 
          }
```

:::{Important}

The paths provided above are just examples. Users should provide `path_to_tracks` and `path_to_he_tracks` based on the location of these folders on their machine.
:::



4. Pass `SSEDict` along with other arguments to the `evolve` function of COSMIC and you are all set. For instance, if you are calculating the evolution of a single binary, the evolve function should look like this: 

```ipython 
bpp, bcm, initC, kick_info = Evolve.evolve(
                              initialbinarytable=single_binary,
                              BSEDict=BSEDict, 
                              SSEDict=SSEDict
                              ) 
``` 

COSMIC will then evolve the system using the provided stellar tracks and store the output in its usual `bpp` and `bcm` arrays.



## BSE 


1. To run METISSE with BSE:

 Get the METISSE-enabled version of BSE through its [GitHub repository](https://github.com/poojanagrawal/BSE-METISSE).


``` console

$ git clone https://github.com/poojanagrawal/BSE-METISSE
$ git submodule update --init --recursive

```
 Alternatively, you can clone BSE and METISSE together: 

``` console

$ git clone --recurse-submodules https://github.com/poojanagrawal/BSE-METISSE.git

```


2. Compile the code using:

``` console

$ ./mk

```
<!-- Make sure use_SSE is false in `bse.f` in the `src` folder before compiling. -->


3. Let us compute the evolution of a single binary using BSE-METISSE, say a 2.5 M<sub>$_\odot$</sub> star with a 1 M<sub>$_\odot$</sub> companion, with an initial orbital period of 100 days in a circular orbit for the evolution time of 13.7 Gyr.

When using BSE, METISSE-specific inputs are read through the `METISSE_input_controls` inlist in the *evolve_metisse.in* file. (See [](input_files.md#metisse-input-controls) for a complete list of input options). 

<!-- **Important:** When using BSE, the file *evolve_metisse.in* should be located in the same directory as the bse executable.  -->


We supply the location of the folders containing [metallicity files](glossary.md#metallicity-file) for hydrogen and helium stars through `METALLICITY_DIR` and `METALLICITY_DIR_HE`. (Make sure to provide relevant paths based on the location of these folders **on your machine** ). 

``` fortran
&METISSE_input_controls

METALLICITY_DIR = '/Users/poojan/Downloads/sample_tracks_solarZ/Hydrogen/'
            
METALLICITY_DIR_HE = '/Users/poojan/Downloads/sample_tracks_solarZ/Helium/'

verbose = .true.

/
 
```

All other inputs including the binary parameters as well as the values of mass and metallicity are read through the BSE input file `binary.in`. 

A typical `binary.in` looks like:


```{code-block} fortran
2.5  1.0 13700 100  1 1 0.02 0.0     ! mass1 mass2 maximum_evolution_time orbital_period_in_days initial_type1 initial_type2 metallicity eccentricity
0.5 0.0 1.0 3.0 0.5                   ! reimers_eta binary_wind_enhancement helium_star_mass_loss_factor alpha_common_envelope lambda_common_envelope
0 1 0 1 0 2 3.0 29769                 ! common_envelope_flag tides_on_off white_dwarf_initial_final_mass_relation white_dwarf_cooling black_hole_kicks neutron_star_remnant_mass max_neutron_star_mass kick_random_seed
0.001 0.01 0.02                        ! timestep_pts1 timesteps_pts2 timestep_pts3
0.0 0.125 1.0 1.5 0.001 10.0 -1.0     ! kick_velocity_dispersion beta_wind_velocity xi_wind_accretion bondi_hoyle_wind_accretion_factor fraction_of_accreted_matter_retained_in_nova_eruption eddington_limit_factor gamma_angular_momentum_loss
1                                      ! sse_flag

```

<!-- ```{code-block}
2.5  1.0 13700 100  1 1 0.02 0.0
0.5 0.0 1.0 3.0 0.5
0 1 0 1 0 2 3.0 29769
0.001 0.01 0.02
0.0 0.125 1.0 1.5 0.001 10.0 -1.0
1
```
where each value in every line is separated by space and corresponds to the following variables.  -->

<!-- 
| Line | Meaning |
|------|---------|
| 1 | `mass1` `mass2` `maximum_evolution_time` `orbital_period_in_days` `initial_type1` `initial_type2` `metallicity` `eccentricity` |
| 2 | `reimers_eta` `binary_wind_enhancement` `helium_star_mass_loss_factor` `alpha_common_envelope` `lambda_common_envelope` |
| 3 | `common_envelope_flag` `tides_flag` `white_dwarf_initial_final_mass_relation_flag` `white_dwarf_cooling_relation` `black_hole_kickflag` `remnant_mass_prescription` `max_neutron_star_mass` `kick_random_seed` |
| 4 | `timestep_pts1` `timesteps_pts2` `timestep_pts3` |
| 5 | `kick_velocity_dispersion` `beta_wind_velocity` `wind_accretion_efficiency` `bondi_hoyle_wind_accretion_factor` `fraction_of_accreted_matter_retained_in_nova_eruption` `eddington_limit_factor_for mass transfer` `gamma_angular_momentum_loss` |
| 6 | `sse_flag: 0 for sse, 1 for metisse` | -->

<!-- :::{Note}
1. If you enter a negative kstar then parameters for an evolved star are required in the order of: current age, initial mass and spin rate, otherwise the star will start on the ZAMS.

2. METISSE currently can not restart the evolution of binaries.

::: -->

4.  Run BSE using:

``` console

$ ./bse_metisse

```

5. BSE will produce the following terminal output:

``` {code-block} console
:class: no-copybutton
METALLICITY_DIR: /Users/poojan/Downloads/sample_tracks_solarZ/Hydrogen

METALLICITY_DIR_HE: /Users/poojan/Downloads/sample_tracks_solarZ/Helium

Input Z is :  2.00000E-02
Reading naked helium star tracks
Found matching Z_files   2.00000E-02
Found           33  tracks.
Minimum initial mass    0.4
Maximum initial mass  148.8
Reading main (hydrogen star) tracks
Found matching Z_files   2.00000E-02
Found           85  tracks.
Minimum initial mass    0.1
Maximum initial mass  299.3
    TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE
    0.0000    2.500    1.000  1  1      137.608  0.00   0.028   0.021  INITIAL 
  567.7662    2.499    1.000  2  1      137.650  0.00   0.064   0.022  KW CHNGE
  574.3388    2.499    1.000  3  1      137.653  0.00   0.120   0.022  KW CHNGE
  580.8607    2.499    1.000  4  1      136.164  0.00   0.524   0.022  KW CHNGE
  773.6237    2.496    1.000  5  1      136.255  0.00   0.295   0.022  KW CHNGE
  782.8053    0.426    1.000 11  1      132.493  0.00   0.000   0.015  KW CHNGE
 9515.7256    0.426    1.000 11  2     1856.760  0.00   0.000   0.002  KW CHNGE
11250.9990    0.426    0.999 11  3     1857.020  0.00   0.000   0.002  KW CHNGE
12178.6533    0.426    0.947 11  4     1927.083  0.00   0.000   0.149  KW CHNGE
12308.4219    0.426    0.935 11  5     1942.590  0.00   0.000   0.023  KW CHNGE
12317.3174    0.426    0.914 11  5     1973.257  0.00   0.000   0.059  BEG SYMB
12317.3223    0.427    0.425 11 11     2795.724  0.00   0.000   0.000  KW CHNGE
13700.0000    0.427    0.425 11 11     3755.136  0.00   0.000   0.000  MAX TIME
```

For more detailed output, check the output file in the `output_binary` folder.

Similar to the standalone mode, to `z_accuracy_limit` can be adjusted in `METISSE_input_controls` to use nearby metallicity value if the required value is not present in the grid. See [](usage_standalone.md#if-the-metallicity-value-is-not-present-in-the-input-grid) for details. 

 <!-- `popbin.f`.  -->


