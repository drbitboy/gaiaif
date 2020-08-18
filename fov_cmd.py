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
  --obspos=X,Y,Z (Observer position, km, triggers parallax correction)
  --obsvel=VX,VY,VZ (Observer velocity, km/s, triggers stellar aberration correction)
  --obsy=YYYY.ddd (Observer time, fractional year, triggers proper mostion correction)

"""

import os
import sys
import pprint
import sqlite3 as sl3
import spiceypy as sp
import gaiaif_util as gifu

dpr = sp.dpr()

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
  obs_pos,obs_vel,obs_year,obs_year_arg = [None]*4
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

    if arg.startswith('--obspos='):
      obs_pos = list(map(float,arg[9:].split(',')))
      continue

    if arg.startswith('--obsvel='):
      obs_vel = list(map(float,arg[9:].split(',')))
      continue

    if arg.startswith('--obsy='):
      obs_year_arg = arg[7:]
      try:
        obs_year = float(obs_year_arg) - 2015.5
      except:
        try:
          obs_year_s,msg = sp.tparse(arg[7:],99)
          assert not msg,'Problem parsing obs_year[{0}]: [{1}]'.format(arg,msg)
        except:
          if '--obsy=' == arg: continue
          raise
        s_2015_5 = '2015-07-02T12:00:00'
        obs_2015_5_s,msg = sp.tparse(s_2015_5,99)
        assert not msg,'Problem parsing obs_year[{0}]: [{1}]'.format(s_2015_5,msg)
        obs_year = (obs_year_s - obs_2015_5_s) / 31557600.0
      continue

    vertex = arg.split(',')
    if 1==len(vertex): vertex = vertex.pop()
    fov_vertices.append(vertex)

    ### End of argument loop

  ### Build build FOV
  fov = gifu.FOV(fov_vertices
                ,obs_pos=obs_pos
                ,obs_vel=obs_vel
                ,obs_year=obs_year
                )

  ### Will need gaialight table if either magnitudes were requested, or
  ### if either parallax or proper motion corrections were requested
  get_light = get_mags or not ((None,None,) == (obs_pos,obs_year,))
  gaiasqls = [GAIASQL(gaia_sl3,mag_min,mag_max,mag_type
                     ,get_ppm,get_light,get_heavy
                     ,*radeclims
                     )
              for radeclims in fov.get_radec_boxes()
             ]

  if do_debug:
    pprint.pprint(locals(),stream=sys.stderr)
    sys.stderr.write('========\n')
    for gaiasql in gaiasqls: sys.stderr.write(gaiasql.query)
    sys.stderr.write('\n========\n')

  rtn_stars = list()

  while len(rtn_stars) < rtn_limit:

    minmag_row = dict(parallax_maspau=None
                     ,pmra_maspy=None
                     ,pmdec_maspy=None
                     )
    minmag_gsql = None

    for row,gsql in [gaiasql.get_row() for gaiasql in gaiasqls]:
      if None is row: continue
      if not (None is minmag_gsql):
        if row['mean_mag'] >= minmag_row['mean_mag']: continue
      minmag_row.update(row)
      minmag_gsql = gsql

    if None is minmag_gsql: break

    star_in_fov,uvstar = fov.star_in_fov([minmag_row['ra'],minmag_row['dec']]
                                        ,parallax_maspau=minmag_row['parallax']
                                        ,pmra_maspy=minmag_row['pmra']
                                        ,pmdec_maspy=minmag_row['pmdec']
                                        )

    if star_in_fov:

      minmag_row['uvstar_corrected'] = list(uvstar)
      (minmag_row['rastar_corrected']
      ,minmag_row['decstar_corrected']
      ,) = sp.vsclg(dpr,sp.recrad(uvstar)[1:3],2)

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
                         ,obs_pos=obs_pos
                         ,obs_vel=obs_vel
                         ,obs_year=obs_year
                         ,obs_year_arg=obs_year_arg
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
    ppm_columns = self.get_ppm and """
      ,gaialight.parallax
      ,gaialight.pmra
      ,gaialight.pmdec""" or ''

    mags_columns = self.get_mags and """
      ,gaialight.phot_g_mean_mag
      ,gaialight.phot_bp_mean_mag
      ,gaialight.phot_rp_mean_mag""" or ''

    heavy_join = self.get_heavy and """
INNER JOIN dbheavy.gaiaheavy
 ON gaiartree.idoffset=dbheavy.gaiaheavy.idoffset
""" or ''

    heavy_columns = self.get_heavy and """
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
      ,gaiaheavy.pmra_pmdec_corr""" or ''

    self.query0 = """
SELECT gaialight.phot_{0}_mean_mag as mean_mag
      ,gaialight.ra as ra
      ,gaialight.dec as dec
      ,gaiartree.idoffset{3}{4}{5}

FROM gaiartree

INNER JOIN gaialight
 ON gaiartree.idoffset=gaialight.idoffset
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
          ,ppm_columns
          ,mags_columns
          ,heavy_columns
          ,heavy_join
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
    if None is self.row:
      rtn = None
    else:
      rtn = dict(zip(self.column_names,tuple(self.row)))
      if 'source_id' in rtn: rtn['source_id'] = str(rtn['source_id'])

    return rtn,self

  def cursor_next(self):
    """Get next .row from cursor"""
    assert not self.done,'Incorrect use of GAIASQL class; contact programmer, code WSNBATGH-GAIASQL-0'
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
