"""fov_cmd.py - command0-line utility to otuput Gaia stars in an FOV

Usage (BASH):

  python fov_cmd.py RA0,Dec0 RA1,Dec1 X2,Y2,Z2

  - Polygonal FOV
  - RA,Dec values in degrees; 0<=RA<360; -90<=Dec<=+90.
  - X2,Y2,Z2 units are not important as long as they are consistent
  - Reference frame assumed to be inertial ICRS; use --j2000 for J2000


  python fov_cmd.py RA0,Dec0 halfang

  - Conical FOV
  - RA0,Dec0 is cone axis
  - halfang is half-angle of cone = angle between axis and cone surface


  python fov_cmd.py RA0,Dec0 RA1,Dec1

  - RA,Dec box
  - RA,Dec values in degrees; 0<=RA<360; -90<=Dec<=+90.
  - N.B. not a real geometric FOV


Output (sys.stdout):

  JSON formatted array, one star per element, by decreasing magnitude

Other options:

  --limit=200
  --mag-max=MaximumMagnitude
  --mag-min=MinimumMagnitude
  --mag-type=g|bp|rp
  --gaia-sqlite3=gaia.sqlite3 (default; also assumes gaia_heavy.sqlite3)
  --j2000 (FOV specified in inertial J2000 instead of ICRS)
  --buffer=pad (around convex FOV or RA,Dec box ;not yet implemented)
  --ppm (Retrieve Parallax and Proper Motions)
  --mags (Retrieve all magnitudes, phot_*_mean_mag)
  --heavy (Retrieve source_id, errors and corr. coeffs, , *error, *corr)

"""

import os
import sys
import pprint
import sqlite3 as sl3
import gaiaif_util as gifu

### Allowed magnitude types
mag_types = set('g bp rp'.split())

do_debug = 'DEBUG' in os.environ


########################################################################
def do_main(argv):

  ######################################################################
  ### Process command-line

  ### - Default values

  rtn_limit = 200
  mag_min,mag_max = None,None
  mag_type = 'g'
  radec_buffer = 0.0   ### Not yet implemented
  gaia_sl3= 'gaia.sqlite3'
  use_j2000 = False
  get_ppm = False
  get_mags = False
  get_heavy = False
  fov_vertices = []

  for arg in argv:
    ### - Loop over arguments

    if arg.startswith('--limit='):
      rtn_limit = int(arg[8:])
      continue

    if arg.startswith('--mag-max='):
      mag_max = float(arg[10:])
      continue

    if arg.startswith('--mag-min='):
      mag_min = float(arg[10:])
      continue

    if arg.startswith('--mag-type='):
      mag_type = arg[11:].strip()
      assert mag_type in mag_types,'Magnitude type argument [{0}] does not specify on of the set of allowed types {1}'.format(arg,mag_types)
      continue

    if arg.startswith('--gaia-sqlite3='):
      gaia_sl3 = arg[15:]
      assert gaia_sl3.endswith('.sqlite3'),'Gaia SQLite3 filepath argument [{0}] does not end in .sqlite3'.format(arg)
      continue

    if '--j2000' == arg:
      use_j2000 = True
      continue

    if '--ppm' == arg:
      get_ppm = True
      continue

    if '--mags' == arg:
      get_mags = True
      continue

    if '--heavy' == arg:
      get_heavy = True
      continue

    if arg.startswith('--buffer='):
      radec_buffer = float(arg[9:])
      sys.stderr.write('Warning:  RA,DEC buffer not yet implemented\n')
      continue

    vertex = arg.split(',')
    if 1==len(vertex): vertex = vertex.pop()
    fov_vertices.append(vertex)

    ### End of argument loop

  ### Build Gaia heavy database path, build FOV
  fov = gifu.FOV(fov_vertices)

  gaiasqls = [GAIASQL(gaia_sl3,mag_min,mag_max,mag_type
                     ,get_ppm,get_mags,get_heavy
                     ,*radeclims
                     )
              for radeclims in fov.get_radec_boxes()
             ]

  if do_debug:
    pprint.pprint(locals(),stream=sys.stderr)
    print('========',file=sys.stderr)
    for gaiasql in gaiasqls: print(gaiasql.query,file=sys.stderr)
    print('========',file=sys.stderr)

  rtn_stars = list()

  while len(rtn_stars) < rtn_limit:

    minmag_row,minmag_gsql = None,None

    for row,gsql in [gaiasql.get_row() for gaiasql in gaiasqls]:
      if None is row: continue
      if not (None is minmag_row):
        if row['mean_mag'] >= minmag_row['mean_mag']: continue
      minmag_row,minmag_gsql = row,gsql

    if None is minmag_row: break

    if fov.star_in_fov([minmag_row['ra'],minmag_row['dec']]):
      rtn_stars.append(minmag_row)

    minmag_gsql.cursor_next()

  return dict(config=dict(limit=rtn_limit
                         ,mag_min=mag_min
                         ,mag_max=mag_max
                         ,mag_type=mag_type
                         ,radec_buffer=radec_buffer
                         ,gaia_sl3=gaia_sl3
                         ,use_j2000=use_j2000
                         ,get_ppm=get_ppm
                         ,get_mags=get_mags
                         ,get_heavy=get_heavy
                         ,fov_vertices=fov_vertices
                         ,fov_type=fov.fovtype
                         )
             ,stars=rtn_stars
             )




########################################################################
class GAIASQL(object):
  def __init__(self,gaia_sl3,lomag,himag,mag_type
              ,get_ppm,get_mags,get_heavy
              ,ralo,rahi,declo,dechi
              ):
    self.gaia_sl3 = gaia_sl3
    self.gaia_heavy_sl3 = '{0}_heavy.sqlite3'.format(gaia_sl3[:-8])
    (self.lomag,self.himag,self.mag_type
    ,self.get_ppm,self.get_mags,self.get_heavy
    ,self.ralo,self.rahi,self.declo,self.dechi
    ,) = (lomag,himag,mag_type
         ,get_ppm,get_mags,get_heavy
         ,ralo,rahi,declo,dechi
         ,)
    assert self.mag_type in mag_types

    self.query_parameters = dict(lomag=self.lomag
                                ,himag=self.himag
                                ,ralo=self.ralo
                                ,rahi=self.rahi
                                ,declo=self.declo
                                ,dechi=self.dechi
                                )
    ppm_columns = """
      ,gaialight.parallax
      ,gaialight.pmra
      ,gaialight.pmdec"""

    mags_columns = """
      ,gaialight.phot_g_mean_mag
      ,gaialight.phot_bp_mean_mag
      ,gaialight.phot_rp_mean_mag"""

    heavy_join = """
INNER JOIN dbheavy.gaiaheavy
 ON gaiartree.offset=dbheavy.gaiaheavy.offset
"""
    heavy_columns = """
      ,gaiaheavy.source_id
      ,gaiaheavy.ra_error
      ,gaiaheavy.dec_error
      ,gaiaheavy.parallax_error
      ,gaiaheavy.pmra_error
      ,gaiaheavy.pmdec_error
      ,gaiaheavy.ra_dec_corr
      ,gaiaheavy.ra_parallax_corr
      ,gaiaheavy.ra_pmra_corr
      ,gaiaheavy.ra_pmdec_corr
      ,gaiaheavy.dec_parallax_corr
      ,gaiaheavy.dec_pmra_corr
      ,gaiaheavy.dec_pmdec_corr
      ,gaiaheavy.parallax_pmra_corr
      ,gaiaheavy.parallax_pmdec_corr
      ,gaiaheavy.pmra_pmdec_corr"""

    self.query0 = """
SELECT gaialight.phot_{0}_mean_mag as mean_mag
      ,gaiartree.ralo as ra
      ,gaiartree.declo as dec
      ,gaiartree.offset{3}{4}{5}

FROM gaiartree

INNER JOIN gaialight
 ON gaiartree.offset=gaialight.offset
{1}
{6}

WHERE gaiartree.ralo <= :rahi
  AND gaiartree.rahi >= :ralo
  AND gaiartree.declo <= :dechi
  AND gaiartree.dechi >= :declo
{2}

ORDER BY phot_{0}_mean_mag

;
""".format(self.mag_type
          ,'{extra_join_light_on}'
          ,'{extra_mag_limits}'
          ,self.get_ppm and ppm_columns or ''
          ,self.get_mags and mags_columns or ''
          ,self.get_heavy and heavy_columns or ''
          ,self.get_heavy and heavy_join or ''
          )

    if None is self.lomag:
      self.extra_join_light_on = ''
      self.extra_mag_limits = ''
    else:
      self.extra_join_light_on = """AND gaiartree.himag >= :lomag\n"""
      self.extra_mag_limits = """  AND gaialight.phot_{0}_mean_mag >= :lomag\n""".format(self.mag_type)

    if not (None is self.himag):
      self.extra_join_light_on += """AND gaiartree.lomag <= :himag\n"""
      self.extra_mag_limits += """  AND gaialight.phot_{0}_mean_mag <= :himag\n""".format(self.mag_type)

    varself = vars(self)
    self.query = self.query0.format(**vars(self))
    self.cursor = sl3.connect(self.gaia_sl3).cursor()
    if self.get_heavy:
      self.cursor.execute("""ATTACH '{0}' as dbheavy""".format(self.gaia_heavy_sl3))
    self.cursor.execute(self.query,self.query_parameters)
    self.column_names = [descs[0] for descs in self.cursor.description]
    self.use___next = hasattr(self.cursor,'__next__')
    self.done = False
    self.count = 0
    self.cursor_next()

  def get_row(self):
    return ((not (None is self.row))
            and dict(zip(self.column_names,tuple(self.row)))
            or None
           ),self

  def cursor_next(self):
    """Get next .row from cursor"""
    assert not self.done,'Bad use of GAIASQL class; contact programmer, code WSNBATGH-GAIASQL-0'
    try:
      if self.use___next: self.row = self.cursor.__next__()
      else              : self.row = self.cursor.next()
      self.count += 1
    except StopIteration as e:
      self.close()
    except:
      raise

  def close(self):
    """Close DB operations; should be called only at StopIteration"""
    self.row = None
    self.done = True
    try: cn = self.cursor.connection
    except: pass
    try: self.cursor.close()
    except: pass
    try: cn.close()
    except: pass


########################################################################
if "__main__" == __name__:

  try: import simplejsonjson as sj
  except: import json as sj
  sj.dump(do_main(sys.argv[1:]),sys.stdout,indent=2)
  sys.stdout.write('\n')
