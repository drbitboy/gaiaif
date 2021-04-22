"""fov_cmd.py - command-line utility to output Gaia stars in an FOV

Usage (Python; see gaiaif(...) below):

  import fov_cmd
  fov_cmd([[RA0,Dec0,],[RA1,Dec1]])  ### Simple RA,DEC window
  fov_cmd([[RA0,DEC0,],halfang)      ### Conical FOV, .3deg half-angle
  fov_cmd([[RAlo,RAhi],[DEClo,DEChi] ### Simple RA,DEC window II*
  fov_cmd([[RA0,Dec0],[RA1,Dec1]
          ,...,[RAn,Decn])           ### Polygonal FOV

  ### * RAlo may be greater than RAhi, so range goes through 360/0

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


  python fov_cmd.py --RAlo,RAhi --DEClo,DEChi

  - RA,Dec box
  - RA,Dec values in degrees; 0<=RA<360; -90<=Dec<=+90.
  - if RAlo < RAhi, then range is [RAlo:RAhi]
  - if RAlo >= RAhi, then range is [RAlo:360] and [0:RAhi]
  - DEClo < DEChi
  - N.B. not a real geometric FOV


Output (sys.stdout):

  JSON formatted array, one star per element, by decreasing magnitude

Other options:

  --ralohi=RAlo,RAhi          ### N.B. RAn,DECn not allowed
  --declohi=DEClo,DEChi       ### N.B. required if --ralohi=... used
  --limit=200
  --magmax=MaximumMagnitude
  --magmin=MinimumMagnitude
  --magtype=g|bp|rp
  --gaia-sqlite3=gaia.sqlite3 (default; also assumes gaia_heavy.sqlite3)
  --j2000 (FOV specified in inertial J2000 instead of ICRS)
  --buffer=pad (around convex FOV or RA,Dec box ;not yet implemented)
  --ppm (Retrieve Parallax and Proper Motions)
  --mags (Retrieve all magnitudes, phot_*_mean_mag)
  --heavy (Retrieve source_id, errors and corr. coeffs, , *error, *corr)
  --obspos=X,Y,Z (Observer position, km, triggers parallax correction)
  --obsvel=VX,VY,VZ (Observer velocity, km/s, triggers stellar aberration correction)
  --obsy=YYYY.ddd (Observer time, fractional year, triggers proper motion correction)

"""

import os
import sys
import pprint
import sqlite3 as sl3
import spiceypy as sp
import gaiaif_util as gifu

dpr = sp.dpr()

### Allowed magnitude types
magtypes = set('g bp rp'.split())

do_debug = 'DEBUG' in os.environ


########################################################################
def gaiaif(fov_vertices
          ,rtn_limit=200
          ,magmin=None,magmax=None
          ,mag_type='g'
          ,radec_buffer=0.0   ### Not yet implemented
          ,gaia_sl3= 'gaia.sqlite3'
          ,j2000=False
          ,ppm=False,mags=False,heavy=False
          ,obs_pos=None,obs_vel=None,obs_year_arg=None
          ,ralohi=[]
          ,declohi=[]
          ,**kwargs
          ):
  """
GAIA R(Tree SQLite3 interface

Sample usage:

  import fov_cmd
  fov_cmd([[10,-45,],[12,-43]])  ### Simple RA,DEC window
  fov_cmd([[10,-45,],0.3)        ### Conical FOV, .3deg half-angle
  fov_cmd([[12,-43,],[10,-43]
          ,[10,-45,],[12,-45])   ### Polygonal FOV

Keywords:

  ,rtn_limit=200
  ,magmax=MaximumMagnitude
  ,magmin=MinimumMagnitude
  ,mag_type='g'|'bp'|'rp'
  ,gaia_sl3='gaia.sqlite3'   Default; also assumes gaia_heavy.sqlite3
  ,j2000=True                FOV is inertial J2000 instead of ICRS
  ,buffer=pad                Pad around convex FOV;not yet implemented
  ,ppm=True                  Retrieve Parallax and Proper Motions
  ,mags=True                 Retrieve all magnitudes, phot_*_mean_mag
  ,heavy                     Retrieve source_id, errors & corr. coeffs
  ,obspos=[X,Y,Z]            Obs posn, km, triggers parallax correction
  ,obsvel=[VX,VY,VZ]         Obs vel, km/s, triggers stellar aberration
  ,obs_year_arg=YYYY.ddd     Obs time, fractional year, triggers PM
  ,obs_year_arg='YYYY-mm-dd-HH:mm:ss.cc'  Alternate obs time, cal date

  """

  if obs_year_arg:
    ### Parse observer time value, convert to y past GAIA DR2 epoch
    try:
      ### Parse fractional year as float
      obs_year = float(obs_year_arg) - 2015.5
    except:
      ### On exception, parse string as calendar date
      obs_year_s,msg = sp.tparse(obs_year_arg,99)
      assert not msg,'Problem parsing obs_year[{0}]: [{1}]'.format(obs_year_arg,msg)
      ### GAIA DR2 epoch
      s_2015_5 = '2015-07-02T12:00:00'
      obs_2015_5_s,msg = sp.tparse(s_2015_5,99)
      assert not msg,'Problem parsing obs_year[{0}]: [{1}]'.format(s_2015_5,msg)
      ### Convert difference to y
      obs_year = (obs_year_s - obs_2015_5_s) / 31557600.0
  else:
    ### No argument supplied
    obs_year = None

  ### Build FOV
  fov = gifu.FOV(fov_vertices
                ,obs_pos=obs_pos
                ,obs_vel=obs_vel
                ,obs_year=obs_year
                ,ralohi=ralohi
                ,declohi=declohi
                )

  ### Will need gaialight table if either proper motions were requested,
  ### or if either parallax or proper motion corrections were requested
  ppm_final = ppm or not ((None,None,) == (obs_pos,obs_year,))
  gaiasqls = [GAIASQL(gaia_sl3,magmin,magmax,mag_type
                     ,ppm_final,mags,heavy
                     ,*radeclims
                     )
              for radeclims in fov.get_radec_boxes()
             ]

  if do_debug:
    pprint.pprint(locals(),stream=sys.stderr)
    sys.stderr.write('========\n')
    for gaiasql in gaiasqls: sys.stderr.write(gaiasql.query)
    sys.stderr.write('\n========\n')

  ### Initialize list of stars that are in FOV
  rtn_stars = list()

  ### Loop over stars
  while len(rtn_stars) < rtn_limit:

    minmag_row = dict(parallax=None
                     ,pmra=None
                     ,pmdec=None
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
      (minmag_row['rastar_delta']
      ,minmag_row['decstar_delta']
      ,) = (minmag_row['rastar_corrected'] - minmag_row['ra']
           ,minmag_row['decstar_corrected'] - minmag_row['dec']
           ,)

      rtn_stars.append(minmag_row)

    minmag_gsql.cursor_next()

  return dict(config=dict(limit=rtn_limit
                         ,magmin=magmin
                         ,magmax=magmax
                         ,mag_type=mag_type
                         ,radec_buffer=radec_buffer
                         ,gaia_sl3=gaia_sl3
                         ,j2000=j2000
                         ,ppm=ppm
                         ,mags=mags
                         ,heavy=heavy
                         ,obs_pos=obs_pos
                         ,obs_vel=obs_vel
                         ,obs_year=obs_year
                         ,obs_year_arg=obs_year_arg
                         ,fov_vertices=fov_vertices
                         ,fov_type=fov.fovtype
                         ,ralohi=fov.ralohi
                         ,declohi=fov.declohi
                         )
             ,stars=rtn_stars
             )


########################################################################
def do_main(argv):
  """Process command line, the call gaiaif(...)"""

  ######################################################################
  ### Process command-line

  fov_vertices = []

  kwargs = dict()

  for arg in argv:
    ### - Loop over arguments

    if arg.startswith('--ralohi='):
      kwargs['ralohi'] = (ralo,rahi,) = list(map(float,arg[9:].split(',')))
      continue

    if arg.startswith('--declohi='):
      kwargs['declohi'] = (declo,dechi,) = list(map(float,arg[10:].split(',')))
      continue

    if arg.startswith('--limit='):
      kwargs['rtn_limit'] = int(arg[8:])
      continue

    if arg.startswith('--magmax='):
      kwargs['magmax'] = float(arg[9:])
      continue

    if arg.startswith('--magmin='):
      kwargs['magmin'] = float(arg[9:])
      continue

    if arg.startswith('--magtype='):
      kwargs['magtype'] = arg[10:].strip()
      assert kwargs['magtype'] in magtypes,'Magnitude type argument [{0}] does not specify on of the set of allowed types {1}'.format(arg,magtypes)
      continue

    if arg.startswith('--gaia-sqlite3='):
      kwargs['gaia_sl3'] = arg[15:]
      assert gaia_sl3.endswith('.sqlite3'),'Gaia SQLite3 filepath argument [{0}] does not end in .sqlite3'.format(arg)
      continue

    if '--j2000' == arg:
      kwargs['j2000'] = True
      continue

    if '--ppm' == arg:
      kwargs['ppm'] = True
      continue

    if '--mags' == arg:
      kwargs['mags'] = True
      continue

    if '--heavy' == arg:
      kwargs['heavy'] = True
      continue

    if arg.startswith('--buffer='):
      kwargs['radec_buffer'] = float(arg[9:])
      sys.stderr.write('Warning:  RA,DEC buffer not yet implemented\n')
      continue

    if arg.startswith('--obspos='):
      kwargs['obs_pos'] = list(map(float,arg[9:].split(',')))
      continue

    if arg.startswith('--obsvel='):
      kwargs['obs_vel'] = list(map(float,arg[9:].split(',')))
      continue

    if arg.startswith('--obsy='):
      kwargs['obs_year_arg'] = arg[7:]
      continue

    vertex = arg.split(',')
    if 1==len(vertex): vertex = vertex.pop()
    fov_vertices.append(vertex)

    ### End of argument loop

  return gaiaif(fov_vertices,**kwargs)




########################################################################
class GAIASQL(object):
  def __init__(self,gaia_sl3,lomag,himag,magtype
              ,ppm,mags,heavy
              ,ralo,rahi,declo,dechi
              ):
    self.gaia_sl3 = gaia_sl3
    self.gaia_heavy_sl3 = '{0}_heavy.sqlite3'.format(gaia_sl3[:-8])
    (self.lomag,self.himag,self.magtype
    ,self.ppm,self.mags,self.heavy
    ,self.ralo,self.rahi,self.declo,self.dechi
    ,) = (lomag,himag,magtype
         ,ppm,mags,heavy
         ,ralo,rahi,declo,dechi
         ,)
    assert self.magtype in magtypes

    self.query_parameters = dict(lomag=self.lomag
                                ,himag=self.himag
                                ,ralo=self.ralo
                                ,rahi=self.rahi
                                ,declo=self.declo
                                ,dechi=self.dechi
                                )
    ppm_columns = self.ppm and """
      ,gaialight.parallax
      ,gaialight.pmra
      ,gaialight.pmdec""" or ''

    mags_columns = self.mags and """
      ,gaialight.phot_g_mean_mag
      ,gaialight.phot_bp_mean_mag
      ,gaialight.phot_rp_mean_mag""" or ''

    heavy_join = self.heavy and """
INNER JOIN dbheavy.gaiaheavy
 ON gaiartree.idoffset=dbheavy.gaiaheavy.idoffset
""" or ''

    heavy_columns = self.heavy and """
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
""".format(self.magtype
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
      self.extra_mag_limits = """  AND gaialight.phot_{0}_mean_mag >= :lomag\n""".format(self.magtype)

    if not (None is self.himag):
      self.extra_join_light_on += """AND gaiartree.lomag <= :himag\n"""
      self.extra_mag_limits += """  AND gaialight.phot_{0}_mean_mag <= :himag\n""".format(self.magtype)

    self.query = self.query0.format(**vars(self))
    self.cursor = sl3.connect(self.gaia_sl3).cursor()
    if self.heavy:
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
