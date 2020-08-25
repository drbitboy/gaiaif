"""
Compare local SQLite3 data and local code against ESA/Gaia TAB web API
for proper motion correction

"""
import os
import sys
import urllib3
import fov_cmd
try: import simplejson as sj
except: import json as sj

httppool = urllib3.PoolManager()

### Query for ESA Gaia TAP web API:
### - two brightest stars;
### - in 2.9-degree circle around RA,Dec=1,2;
### - proper motion correction to time 2016.5, 1y post Gaia DR2 epoch.

QUERY="""
SELECT source_id
      ,ra
      ,dec
      ,array_element(arr,1) as rastar_corrected
      ,array_element(arr,2) as decstar_corrected
      ,pmra
      ,pmdec
      ,parallax
      ,phot_g_mean_mag as mean_mag

FROM (
  SELECT TOP 2
  source_id
 ,ra
 ,dec
 ,pmra
 ,pmdec
 ,parallax
 ,phot_g_mean_mag
 ,epoch_prop(ra,dec,parallax,pmra,pmdec,radial_velocity,2015.5,2016.5) as arr

  FROM gaiadr2.gaia_source

  WHERE 1=CONTAINS(point('ICRS',ra,dec)
                  ,CIRCLE('ICRS',1,2,2.9)
                  )

  ORDER BY PHOT_G_MEAN_MAG
) as p
;"""

def intorflt(s):
  """Convert string to int or to float"""
  try: return int(s)
  except: return float(s)

class GENERIC(dict):
  """Dict class with attributer for each key"""
  def __init__(self,names,vals):
    self.update(dict(zip(names,map(intorflt,vals))))
    for name in self:
      setattr(self,name,intorflt(self[name]))

def get_tap_data():
  """Query data from ESA Gaia TAP web api; return list of dicts"""
  r = httppool.request('POST'
                      ,"https://gea.esac.esa.int/tap-server/tap/sync"
                      ,fields=dict(FORMAT='JSON'
                                  ,LANG='ADQL'
                                  ,REQUEST='doQuery'
                                  ,QUERY=QUERY
                                  )
                      )

  rtn_dict = sj.loads(r.data.decode('utf-8'))
  names = [d['name'] for d in rtn_dict['metadata']]
  rtn_dict['stars'] = [GENERIC(names,row)
                       for row in rtn_dict['data']
                      ]

  return rtn_dict

if "__main__" == __name__:
  import pprint

  ### Query TAP server
  tap_data = get_tap_data()

  if '--test-gaia' in sys.argv[1:]:

    ### Make similar query local data from SQLite3 database (DB)
    sl3_data = fov_cmd.do_main('1,2 2.9 --mag-type=g --limit=2 --obsy=2016.5'.split())

    ### Get GENERIC class instances from each
    tap_stars,sl3_stars = [d['stars'] for d in (tap_data,sl3_data,)]

    ### Store differences greater than 2e-12 in values of matching keys
    keys = set(tap_stars[0].keys()).intersection(set(sl3_stars[0].keys()))
    diffs=[dict(zip(keys,[sl3row[key]-taprow[key]
                          for key in keys
                          if abs(sl3row[key]-taprow[key])>2e-12
                         ]
                   )
               )
           for taprow,sl3row in zip(tap_stars,sl3_stars)
          ]

    ### Test results:  all rows should be empty dicts
    try:
      for row in diffs: assert not row
    except:
      pprint.pprint(dict(diffs=diffs))
      raise

  else:
    pprint.pprint(tap_data)
