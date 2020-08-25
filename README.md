# Interfaces to Gaia Data Release 2 star catalog SQLite3 R-Tree database

## Cf. https://github.com/drbitboy/Tycho2_SQLite_Rtree/tree/master/gaia

## Sample database files of a very small subset of Gaia catalog are at https://github.com/drbitboy/gaia_subset

    Sample Usage:

      octave:1> addpath('jsonlab')     %%% Cf. git clone https://github.com/fangq/jsonlab.git jsonlab/
      octave:2> x.fov={{1,2},.3};      %%% Configure 0.3deg radius circular FOV at RA=1deg, Dec=2deg
      octave:3> x.limit=2;             %%% Configure to return only two lowest (G) magnitude stars
      octave:4> fov_cmd(x)             %%% Get the data ...

      ans = 
      {
        [1,1] =

          scalar structure containing the fields:

            mean_mag =  7.9053
            ra =  1.1257
            dec =  2.2674
            offset =  116655034

        [1,2] =

          scalar structure containing the fields:

            mean_mag =  9.9577
            ra =  1.0022
            dec =  2.2425
            offset =  116655047

      }

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    See also comments in fov_cmd.m and test_octave target in Makefile

## Manifest

* README.md - this file
* Makefile - make file to run tests
* fov_cmd.m - Matlab/Octave script to test FOV interface
* fov_cmd.py - Python script called by fov_cmd.m
* gaia_icrs_fk.py - Transformation between J2000 and ICRS reference frames
* gaiaif_util.py - Gaia interface utilities; implements FOV and FOVSIDE classes
* test_parallax_stellar_aberr.py - Script to test stellar aberration calculation
* test_query.py - sample query of Gaia SQLite3 database (DB)
* urlget_test.py - Compare local Gaia interface against ESA/Gaia TAB web API
* validate_delta_ra_formula.py - Validate formula that calculates half-RA (Right Ascension) difference of two planes that contain a conical FOV (Field Of View)
* validate_gaiaif_fov.py - Test code for gaiaif_util.FOV class
* gaiaif.py - Initial attempt at Gaia interface; not yet finished; use gaiaif_util.py instead
* 00readme.txt - Text version of Sample usage above
* Extra-repo files
  * Makefile.dotar - Include file for Makefile
  * make_test_octave.log - Output of [make test_octave]
  * gaia.sqlite3 - Typically a symlink to Gaia SQLite3 light DB
  * gaia_heavy.sqlite3 - Typically a symlink to Gaia SQLite3 heavy DB
  * de421.bsp - Solar system ephemeris downloaded and used by Makefile tests
  * jsonlab/ - Location of JSON parsing code for fov_cmd.py; Git-cloned by Makefile
