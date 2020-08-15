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

    See also comments in fov_cmd.m
