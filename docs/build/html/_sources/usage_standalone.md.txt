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
Input to METISSE in the standalone mode is provided through two [Fortran namelists](glossary.md#fortran-namelists): `SSE_input_controls` and `METISSE_input_controls`. 
See sections [](input_files.md#sse-input-controls) and [](input_files.md#metisse-input-controls) for a complete list of input options. 


`SSE_input_controls` is contained in the file called *main.input* and
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
In the standalone mode of METISSE, `METISSE_input_controls` is contained in the file called *metisse.input*. 


We use the variable `METALLICITY_DIR` to supply paths to the folder containing [metallicity files](glossary.md#metallicity-file) for normal hydrogen stars and the variable `METALLICITY_DIR_HE` for naked helium stars. For the pre-packaged grid of stellar tracks available with METISSE, this is the path to the `hydrogen` and `helium` folders respectively.


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
 count        1  input mass =   1.000
     Time        0.0    Phase         MS     Mass   1.000
     Time     9515.7    Phase         HG     Mass   1.000
     Time    11251.0    Phase        FGB     Mass   0.999
     Time    12000.0    Phase        FGB     Mass   0.999
 -------------------------------------------------------------------------
 Reached the end of the program

```

Since we had set, `write_output_to_file = .true.` in `SSE_input_controls`, a [SSE-style](glossary.md#sse-style-file) output file named 'evolve_00100M.dat' will also be generated in the *output* directory, containing a more detailed evolutionary history of the star.


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

We also set `verbose = .false.` in metisse.input to avoid flooding our screen with evolution details of 10,000 stars. 

Running METISSE does not produce any on-screen output this time. 

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
 Input Z is :  1.60000E-02
 Reading naked helium star tracks
 No matching Z_files found with Z_accuracy_limit =  1.00000E-02
 Switching to SSE formulae for helium stars 
 Reading main (hydrogen star) tracks
 No matching Z_files found with Z_accuracy_limit =  1.00000E-02
 If needed, Z_accuracy_limit can be increased to match one of the available Z_files 
  1.00000E-05  1.40000E-05  1.90000E-05  2.00000E-05  2.60000E-05  3.60000E-05  4.90000E-05  6.70000E-05  9.20000E-05  1.00000E-04  1.27000E-04  1.40000E-04  1.74000E-04  2.00000E-04  2.40000E-04  3.29000E-04  4.52000E-04  6.21000E-04  8.53000E-04  1.00000E-03  1.17200E-03  1.40000E-03  1.61000E-03  2.00000E-03  2.21200E-03  3.03900E-03  4.17500E-03  5.73600E-03  7.88000E-03  1.00000E-02  1.08260E-02  1.40000E-02  1.48740E-02  2.00000E-02  2.04340E-02  2.80720E-02  3.85660E-02  5.29830E-02  7.27900E-02
STOP Fatal error: terminating METISSE

````

METISSE does not interpolate in metallicity. So ideally, we should compute input tracks with Z = 0.016 and provide them to METISSE. 

Another option is that we can ask METISSE to use nearby metallicity values. 
This is achieved by increasing `Z_accuracy_limit` in `METISSE_input_controls`.
METISSE lists available metallicity values in the grid, whenever it cannot find the input value in the grid. We can increase `Z_accuracy_limit` to `8d-2` so that METISSE can run using the closest metallicity in the list, which is `0.0148`. 

``` console
$ ./metisse
 Input Z is :  1.60000E-02
 Reading naked helium star tracks
 Found matching Z_files   1.48740E-02
 Found           35  tracks.
 Minimum initial mass    0.4
 Maximum initial mass  148.7
 Reading main (hydrogen star) tracks
 Found matching Z_files   1.48740E-02
 Found           85  tracks.
 Minimum initial mass    0.1
 Maximum initial mass  299.4
 count        1  input mass =   1.000
     Time        0.0    Phase         MS     Mass   1.000
     Time     8024.2    Phase         HG     Mass   1.000
     Time     9563.4    Phase        FGB     Mass   0.999
     Time    10340.5    Phase       CHeB     Mass   0.958
     Time    10461.4    Phase       EAGB     Mass   0.955
     Time    10471.1    Phase      CO_WD     Mass   0.421
     Time    12000.0    Phase      CO_WD     Mass   0.421
 -------------------------------------------------------------------------
 Reached the end of the program

```

METISSE now runs but for metallicity `0.0148`. 
