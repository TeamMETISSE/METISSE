# Using METISSE to evolve single stars

 
METISSE can be used independently to compute the evolution of one star or a population of single stars. 


In order to run METISSE, we first need to compile it. To do so, open a command line shell and inside the METISSE folder execute:

```console

$ ./mk

```

You will now see an executable named `metisse` which was not there before. 

:::{Note}

The code does not need to be re-compiled unless you make changes inside the source `src` directory. 

:::


## Evolving one star 

Let us compute the evolution of a star with initial mass 1 M<sub>$_\odot$</sub> star, metallicity `Z = 0.02` up to the age of 12 Gyr. 
Input to METISSE in the standalone mode is provided through two [Fortran namelists](acronyms_definitions.md#fortran-namelists): `SSE_input_controls` and `METISSE_input_controls`. 
See sections [](input_files.md#sse-input-controls) and [](input_files.md#metisse-input-controls) for a complete list of input options. 


`SSE_input_controls` is contained in the file called `main.input` and
is used to provide evolution details of the star.
 
```
&SSE_input_controls

initial_Z = 0.02

max_age = 1.2d4   

number_of_tracks = 1
min_mass = 1.d0

write_output_to_file = .true.   
/

```

We also need to provide the details of the input tracks to METISSE through the namelist `METISSE_input_controls`. 
In the standalone mode of METISSE, `METISSE_input_controls` is contained in the file called `metisse.input`. 


We use the variable `METALLICITY_DIR` to supply paths to the folder containing [metallicity files](acronyms_definitions.md#metallicity-file) for normal hydrogen stars and the variable `METALLICITY_DIR_HE` for naked helium or stripped stars. For the pre-packaged grid of stellar tracks available with METISSE, this is the path to the `hydrogen` and `helium` folders respectively.


```
&METISSE_input_controls


METALLICITY_DIR = '/Users/poojan/stellar_tracks/MESA/big_z/hydrogen'
            
METALLICITY_DIR_HE = '/Users/poojan/stellar_tracks/MESA/big_z/helium'

verbose = .true.

/

```

:::{Important}

The paths provided above are just examples. Users should input `METALLICITY_DIR` and `METALLICITY_DIR_HE` based on the actual location of these folders on their machine after downloading the grid.
:::


Note that we have also set `verbose` to true so that we can see details of the input tracks on the screen. 
If verbose is False, then these details are saved in the `tracks_log.txt` file in the METISSE directory. 


To run METISSE simply type `./metisse` on the command line and hit enter. METISSE will produce the following output on the screen:

``` console

 $ ./metisse
 Input Z is :  0.020000
 Reading naked helium star tracks
 Reading input files from: /Users/poojan/stellar_tracks/MESA/big_z/helium/./folder_Z0020000/eeps
 Found:           34  tracks
 Minimum initial mass    0.4
 Maximum initial mass  148.8
 Reading main (hydrogen star) tracks
 Reading input files from: /Users/poojan/stellar_tracks/MESA/big_z/hydrogen/./folder_Z0020000/eeps
 Found:           85  tracks
 Minimum initial mass    0.1
 Maximum initial mass  299.3
 count        1  input mass =   1.000
     Time        0.0    Phase         MS     Mass   1.000
     Time     9515.7    Phase         HG     Mass   1.000
     Time    11251.0    Phase        FGB     Mass   0.999
     Time    12000.0    Phase        FGB     Mass   0.999
 -------------------------------------------------------------------------
 Reached end of the program

```

Since we had set, `write_output_to_file = .true.` in `SSE_input_controls`, a [SSE-style](acronyms_definitions.md#sse-style-file) output file named 'evolve_00100M.dat' will also be generated in the *output* directory, containing a more detailed evolutionary history of the star.


## Evolving a stellar population

To compute the evolution of a single stellar population containing, for example, 10,000 stars with initial masses uniformly distributed between 1 and 100 M<sub>$_\odot$</sub>, we modify the `SSE_input_controls` from the single star case as shown below. 


```
&SSE_input_controls

initial_Z = 0.02

max_age = 1.2d4   

number_of_tracks = 1d+4
min_mass = 1.d0

max_mass = 100.d0

write_track_to_file = .true.   
/

```

We also set `verbose = .false.` in  `metisse.input` to avoid flooding our screen with evolution details of 10,000 stars. 

Running METISSE does not produces any on-screen output this time. 

``` console

$ ./metisse

```

We also kept `write_track_to_file = .true.` in `SSE_input_controls`. Therefore, the `output` directory will now contain the evolutionary histories of all 10,000 stars. 


To compute a population of single stars with any other distribution of initial masses, we can list initial masses in a file provide masses in a text file (one per line) and provide the location of that file through `input_mass_file`. 
If the masses are listed in a file called `my_custom_distribution.txt`, then `SSE_input_controls` will look like this:


```
&SSE_input_controls

initial_Z = 0.02

max_age = 1.2d4   

number_of_tracks = 1d+4

read_mass_from_file = .true.  

input_mass_file = 'my_custom_distribution.txt'

write_track_to_file = .true.  

/

```

Note that we still need to define `number_of_tracks` as METISSE will read `input_mass_file` for that many stars and will raise errors if the file contains fewer stars.  

## If the metallicity value is not present in the input grid


Let's say in the above example of a single star, we set `initial_Z = 0.016`. 
Running METISSE now gives an error, as the metallicity value is not present in the grid of input tracks.

``` console
$ ./metisse 
 Input Z is :  0.016000
 Reading naked helium star tracks
 METISSE error: metallicity value =   1.6000000000000000E-002 not found amongst given Z_files
 Check metallicity_file_list and value of Z_files for each file
 If needed, Z_accuracy_limit can be relaxed (set to a greater value).
 Switching to SSE formulae for helium stars 
 Reading main (hydrogen star) tracks
 METISSE error: metallicity value =   1.6000000000000000E-002 not found amongst given Z_files
 Check metallicity_file_list and value of Z_files for each file
 If needed, Z_accuracy_limit can be relaxed (set to a greater value).
STOP Fatal error: terminating METISSE

````

METISSE does not interpolate in metallicity. So ideally, we should compute input tracks with Z = 0.016 and supply them to METISSE. 

Another option is that we can ask METISSE to use nearby metallicity values. 
This is achieved by the variable `Z_accuracy_limit` in `METISSE_input_controls`. Be default it is set to `1d-2`. We increase it to `8d-2`. 

``` console
$ ./metisse
 Input Z is :  0.016000
 Reading naked helium star tracks
 Reading input files from: /Users/poojan/stellar_tracks/MESA/big_z/helium/./folder_Z0010000/eeps
 Found:           37  tracks
 Minimum initial mass    0.4
 Maximum initial mass  126.7
 Reading main (hydrogen star) tracks
 Reading input files from: /Users/poojan/stellar_tracks/MESA/big_z/hydrogen/./folder_Z0010000/eeps
 Found:           84  tracks
 Minimum initial mass    0.1
 Maximum initial mass  299.7
 count        1  input mass =   1.000
     Time        0.0    Phase         MS     Mass   1.000
     Time     6442.1    Phase         HG     Mass   1.000
     Time     7754.5    Phase        FGB     Mass   0.999
     Time     8393.0    Phase       CHeB     Mass   0.961
     Time     8513.5    Phase       EAGB     Mass   0.958
     Time     8521.9    Phase      CO_WD     Mass   0.432
     Time    12000.0    Phase      CO_WD     Mass   0.432
 -------------------------------------------------------------------------
 Reached end of the program

```

METISSE now runs using the closest metallicity value it can find. In this case, it is `0.01`. 
