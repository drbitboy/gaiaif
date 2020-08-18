"""
gaiaif.py - Interface to Gaia data in SQLite3 database files

"""
import re
import os
import sys
import spiceypy as sp

try:
  assert ALLRTREE
except:
  ###import traceback
  ###traceback.print_exc()
  ALLRTREE,ALLLIGHT,ALLHEAVY = all_table_keys = """lgt.gaiartree lgt.gaialight hvy.gaiaheavy""".strip().split()
  table_columns = { ALLRTREE: """ralo,rahi,declo,dechi,lomag,himag""".strip().split(',')
                  , ALLLIGHT: """ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag""".strip().split(',')
                  , ALLHEAVY: """source_id,ra_error,dec_error,parallax_error,pmra_error,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr""".strip().split(',')
                  }

  for table_key in all_table_keys:
    table_columns.update(dict([(column.strip(),table_key,) for column in table_columns[table_key]]))

  class GaiaUseException(Exception): pass

  CIRCLETYPE,RADECBOXTYPE,POLYGONTYPE = 'circle radecbox polygon'.split()

  J2000,ICRS = 'J2000 ICRS'.split()
  Reffrms = dict()
  for value in (J2000,ICRS,):
    for key in (value,value.lower(),value[0],value[0].lower(),):
      Reffrms[key] = value


def generic_rdm_query(return_columns
                     ,where_limits
                     ,count_limit
                     ):
  """
Generic query within RA,Dec,Magnitude limits

"""


def fov_query(*args):
  """
User query for Gaia stars in FOV

Usage:

  Descr,data = fov_query(ppfx,rfrm,fov[,mags[,cols[,ct[,epch[,opos[,ovel]]]]])

  Return

    Descr:   Names of columns of data
    values:  2-D array of returned rows of data


  Arguments; see Input Argument Details below for more information

    ppfx:  Path prefix to SQLite3 files, without [_heavy].sqlite3 suffix
    rfrm:  reference frame of user inputs and return data, J2000 or ICRS
    fov:   FOV [RA,Dec] values to constrain lookup
    mags:  Magnitude constraint(s)
    cols:  which columns' values to return
    ct:    Maximum count of stars to return, plus sorting specification
    epch:  Epoch for proper motion; None if not desired
    opos:  Observer position, for parallax correction
    ovel:  Observer velocity, for stellar aberration correction

    - Any arguments after FOV are optional, but all, before the last one
      present, must be supplied for that last one to be parsed correctly


  Input Arguments Details

    ppfx:  String, or pair of strings
           Single string to SQLite3 light database (DB) file; heavy file
           will be build by inserting _heavy before the last full-stop
           e.g.
           - 'sub-directory/gaia.sqlite3'
             - Light DB will be  sub-directory/gaia.sqlite3
             - Heavy DB will be  sub-directory/gaia_heavy.sqlite3
               - Split name before last full-stop

           Pair of strings to SQLite3 database (DB) files e.g.
           - ['sub-dir/gaia.sqlite3','sub-dir/gaiahvy.db']
             - Light DB will be  sub-directory/gaia.sqlite3
             - Heavy DB will be  sub-directory/gaiahvy.db

           N.B. if the heavy DB is not required by the chosen input
                arguments, the name will be built but it will not be
                opened, so the heavy file does not have to exist.

    rfrm:  String

           The Solar System Barycenteric reference frame of the inputs
           and the returned outputs (RA, Dec, position, velocity, etc.)
           Will be either 'J2000' or 'ICRS'; acceptable alternate
           shorthand:  gaiaif.J2000 'J2000' 'j2000' 'J' 'j'
                       gaiaif.ICRS 'ICRS' 'icrs' 'I' 'i'

    fov:   List
           Specification of FOV shape to spatailly constrain lookup
           - [[RA,Dec],hang] cone with axis=[RA,Dec] and half-angle=hang
           - [[RA,Dec],[RA,Dec]] RA,Dec box opposite corners
           - [[RA,Dec],[RA,Dec],...,[RA,Dec]] RA,DEC polygon vertices

    mags:  List, or None
           Magnitude constraint(s):  [upper mag]; [upper mag, lower mag]
           - or use None for no magnitude constraint

    cols:  List, or None
           Names of first columns of data to return N.B. may be extended
           e.g. parallax correction.  Possible values:

           From R-Tree table (does not require gaia_heavy.sqlite3):
           - ralo or rahi:  nominal RA at Gaia epoch (2015.5), deg
           - declo or dechi:  nominal Dec at Gaia epoch (2015.5), deg
           - lomag:  low magnitude
           - himag:  high magnitude
           - gaiaif.ALLRTREE:  all of the R-Tree table columns

           From Light table (does not require gaia_heavy.sqlite3):
           - parallax:  absolute stellar parallax (mas)
           - pmra:  proper motion in RA (mas/y)
           - pmdec:  proper motion in Dec (mas/y)
           - phot_g_mean_mag:  mean magnitude, G band
           - phot_bp_mean_mag:  mean magnitude, BP band
           - gaiaif.ALLLIGHT:  all of the Light table columns

           From Heavy table (requires gaia_heave.sqlite3):
           - source_id:  Gaia source ID, 64-bit integer
           - ra_error:  RA error, mas
           - dec_error:  Dec error, mas
           - parallax_error:  Dec error, mas
           - pmra_error:  RA proper motion error, mas/y
           - pmdec_error:  Dec proper motion error, mas/y
           - X_Y_corr:  Corr. coeff btw X and Y, dimensionless, [-1,+1]
             - X is one of ra, dec, parallax, and pmra
             - Y is dec, parallax, pmra and  pmdec
             - X is not the same as Y
           - gaiaif.ALLHEAVY:  all of the Heavy table columns

           If None, defaults to gaiaif.RTREE (RA, Dec, lomag, himag)

    ct:    List, or integer, or None
           Maximum count of stars to return and sorting specification

           [number, 'column sort-order']
           - Sort-order is either ASC or DESC, case insensitive

           number
           - Limit with no sorting, no limit

           None
           - no sorting, no limit

           e.g.
           - [100,'lomag asc']
           - [200,'phot_g_mean_mag asc']
           - [300,'declo desc']
           - 400
           
    epch:  String, or float, or None
           Epoch for proper motion correction

           UTC, String
           - e.g. '2021-01-23T34:56:01.23456'
           - Requires SPICE LEAPSECOND Kernel to be FURNSHed (LSK)

           TDB, float, integer or other numeric
           - s past J2000 epoch, float

           None
           - if no proper motion correction is desired

    opos:  3-Vector sequence, or numpy array, or None
           Observer position, for parallax correction, km

           If None, no parallax correction

    ovel:  3-Vector sequence or numpy array
           Observer velocity, for stellar aberration correction, km/s

           If None, no aberration correction

"""

  ### Extract and parse input arguments

  ### 0 ppfx:  nothing to do; this will be passed to the DB call

  ### 1 rfrm:  reference frame
  rfrm = Reffrms[args[1]]

  ### 2 fov:  [[RA,Dec],hang] or [[RA,Dec],[RA,DEC]...]
  fovvertices = args[2]
  assert len(fovvertices)>1
  ra0,dec0 = map(float,fovvertices[0])
  assert ra0<=360.0 and ra0>=0.0
  assert dec0<=90.0 and dec0>=-90.0
  try:
    hang = float(fovvertices[1])
    fovtype = CIRCLETYPE
  except:
    for pair in fovvertices[1:]:
      ra,dec = map(float,pair)
      assert ra<=360.0 and ra>=0.0
      assert dec<=90.0 and dec>=-90.0
    fovtype = 2==len(fovvertices) and RADECBOXTYPE or POLYGONTYPE


if "__main__"==__name__:
  print(table_columns)
  generic_rdm_query(None,None,None)

