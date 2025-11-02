
# METISSE with COSMIC 

COSMIC [(Breivik et al. 2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract) is a binary population synthesis code designed to model the evolution of stellar systems across a wide range of astrophysical environments. It builds upon the rapid binary evolution algorithms of BSE, extending them with population synthesis capabilities and easy to use Python framework. 
METISSE can be used as the stellar evolution module inside COSMIC.

This section provides a short, step-by-step guide to using METISSE with COSMIC.
For more details please refer to COSMIC's [official documentation](https://cosmic-popsynth.github.io/COSMIC/index.html). 

1. Follow COSMIC’s official [installation guide](https://cosmic-popsynth.github.io/COSMIC/pages/install.html) to load the COSMIC environment and activate it in your terminal.

    If you don’t already have COSMIC installed, download it from [GitHub](https://github.com/COSMIC-PopSynth/COSMIC/tree/develop).


    ``` console

    $ git clone https://github.com/COSMIC-PopSynth/COSMIC

    ```

    Next, change into the COSMIC directory and compile the library: 


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

    METISSE is now ready to use with COSMIC.

<!-- 1. Get COSMIC (version 4.0)

    Follow COSMIC’s official [installation guide](https://cosmic-popsynth.github.io/COSMIC/pages/install.html) to load the COSMIC environment and activate it in your terminal.

    If you don’t already have COSMIC installed, download it from [GitHub](https://github.com/COSMIC-PopSynth/COSMIC/tree/develop).

    ``` console

    $ git clone --recurse-submodules https://github.com/COSMIC-PopSynth/COSMIC

    ```

    `--recurse-submodules` downloads METISSE along with COSMIC.

    If you are an existing COSMIC user, update your repository and fetch the METISSE submodule:

    ``` console

    $ git submodule update --init --recursive

    ```

2. Next, change into the COSMIC directory and compile the library: 


    ``` console

    $ pip install .

    ```

    METISSE is now ready to use with COSMIC.  -->


3. Configure METISSE inside COSMIC

    We will try evolving an example binary following the instructions on [this page](https://cosmic-popsynth.github.io/COSMIC/pages/examples.html). 

    The usual input namelists for METISSE are not read when using COSMIC. Instead a Python dictionary, `SSEDict` is used to set the `stellar_engine` and provide the location of the folders containing {term}`Metallicity File` for hydrogen and helium stars through variables `path_to_tracks` and `path_to_he_tracks`.

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

4. The various binary evolution prescriptions used by BSE are set using a Python dictionary, `BSEDict`:

    ``` ipython

    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': [1.0,1.0], 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : [-1,-1], 'rtmsflag' : 0, 'wd_mass_lim': 1}

    ```

5. You can define the binary’s initial properties using the InitialBinaryTable class:

    ``` ipython
    single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757,
                                                    ecc=0.448872, tphysf=13700.0,
                                                    kstar1=1, kstar2=1, metallicity=0.002)
    ```

6. Finally, evolve the binary by passing `SSEDict`, `BSEDict` and `initialbinarytable` to the `evolve` function of COSMIC.

    ```ipython 
    bpp, bcm, initC, kick_info = Evolve.evolve(
                                initialbinarytable=single_binary,
                                BSEDict=BSEDict, 
                                SSEDict=SSEDict
                                ) 
    ``` 

    COSMIC will then evolve the system using the provided stellar tracks and store results in its standard output arrays:

    `bpp` → binary properties at important stages in the binary’s evolution 

    `bcm` → binary parameters at user specified time steps during the binary’s evolution.

