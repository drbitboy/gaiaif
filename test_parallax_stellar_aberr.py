"""
Compare gaiaif_util.py corrections for parallax and stellar aberration
to equivalent correction in SPICE

"""
import os
import math
import fov_cmd
import sqlite3 as sl3
import spiceypy as sp
import gaiaif_util as gifu

### Setup debugging and conversions; define some string constants
do_debug = 'DEBUG' in os.environ
rpd,aupkm = sp.rpd(),sp.convrt(1.,'km','au')
(earth,ssb,j2000,none,LT,LTS,ra_dec,declination,equals
,) = 'earth 0 j2000 none lt lt+s ra/dec declination ='.upper().split()

sp.furnsh('de421.bsp')

gaia_db = 'gaia.sqlite3'

### Find ET of spring equinox wrt Solar System Barycenter (SSB) & earth
cnfine=sp.utils.support_types.SPICEDOUBLE_CELL(2)
result=sp.utils.support_types.SPICEDOUBLE_CELL(200)
sp.wninsd(0.,sp.pi()*.5e7,cnfine)
r=sp.gfposc(ssb,j2000,none,earth,ra_dec,declination,equals,0.,0.,10*sp.spd(),200,cnfine,result)

### Get SSB->Earth and Earth->SSB vectors at that time
et = sp.wnfetd(result,0)[0]
ssb2e = sp.spkezr(earth,et,j2000,none,ssb)[0]
e2ssblt = sp.spkezr(ssb,et,j2000,LT,earth)[0]
if do_debug:
  ### N.B. Light-time correction on Earth->SSB will have no effect
  print(dict(etcal=sp.etcal(et,99)
            ,ssb2earth_au=sp.vscl(aupkm,ssb2e[:3])
            ,earth2ssb_au=sp.vscl(aupkm,e2ssblt[:3])
            ,))
  assert 0.0==sp.vnormg(sp.vaddg(ssb2e,e2ssblt,6),6)

### Get a star (ra,dec,parallax) near the pole with parallax > 10mas/y
cn=sl3.connect(gaia_db)
cu=cn.cursor()
cu.execute("""
SELECT gr.idoffset
      ,gl.ra
      ,gl.dec
      ,gl.parallax
FROM gaiartree as gr
INNER JOIN gaialight as gl
ON gr.idoffset=gl.idoffset
WHERE gr.dechi> 89.5
  AND gl.parallax>10.0
ORDER BY gl.dec DESC
LIMIT 1
;""")

idoffset,ra,dec,parallax = result = list(cu.fetchone())

if do_debug: print(dict(zip([d[0] for d in cu.description],result)))

### Setup one gaiaif_util.FOV instance with only a parallax correction,
### and another with parallax and stellar aberration corrections
### - The FOV is large here as the objective is the correction, not the
###   True/False result of FOV.star_in_fov
fov_noab = gifu.FOV([[0.,90.],5.],obs_pos=ssb2e[:3])
fov_ab = gifu.FOV([[0.,90.],5.],obs_pos=ssb2e[:3],obs_vel=ssb2e[3:])

### Call .star_in_fov with the current near-pole star to get corrected
### positions
in_noab,uvstar_noab = fov_noab.star_in_fov([ra,dec],parallax_maspau=parallax)
in_ab,uvstar_ab = fov_ab.star_in_fov([ra,dec],parallax_maspau=parallax)

### Ensure both stars were in the FOVs
assert in_noab and in_ab

### Calculate range to the star from the parallax value
range2star = 1.0/(aupkm * math.tan(parallax*rpd/3.6e6))
vstar = sp.radrec(range2star,ra*rpd,dec*rpd)

### Use SPICE to calculate vector from earth to the the SSB-relative
### fixed star position without any correction
earth2star_none = sp.spkcpt(vstar,ssb,j2000
                           ,et,j2000,"observer",none
                           ,"earth"
                           )[0]

### Use SPICE to calculate vector from earth to the the SSB-relative
### fixed star position with LT+S correction
earth2star_lts = sp.spkcpt(vstar,ssb,j2000
                          ,et,j2000,"observer",LTS
                          ,"earth"
                          )[0]

if do_debug:
  ### Ensure that light-time correction is zero for earth->star vector
  ### for star at fixed positon relative to SSB
  earth2star_lt   = sp.spkcpt(vstar,ssb,j2000
                             ,et,j2000,"observer",LT
                             ,"earth"
                             )[0]
  assert 0.0 == sp.vnormg(sp.vsubg(earth2star_none
                                  ,earth2star_lt
                                  ,6),6)

### Calculate differences between gaiif_util.py and SPICE results ...
ab_vec_spice = sp.vsub(earth2star_lts[:3],earth2star_none[:3])
ab_vec_gaia = sp.vscl(range2star,sp.vsub(uvstar_ab,uvstar_noab))
ab_vec_diff = sp.vsub(ab_vec_spice,ab_vec_gaia)
ab_diff_frac = sp.vnorm(ab_vec_diff) / sp.vnorm(ab_vec_spice)

no_ab_diff = sp.vsep(uvstar_noab,sp.vhat(earth2star_none[:3]))
ab_diff = sp.vsep(uvstar_ab,sp.vhat(earth2star_lts[:3]))
ab_diff_frac = sp.vsep(uvstar_ab,sp.vhat(earth2star_lts[:3]))

try:
  ### ... Success
  assert 1e-15 > no_ab_diff
  assert 1e-8 > ab_diff
except:
  ### ... Failure
  print(dict(no_ab_diff=no_ab_diff
            ,ab_diff=ab_diff
            ,ab_diff_frac=ab_diff_frac
            ,ab_mag_spice=sp.vsep(earth2star_none[:3],earth2star_lts[:3])
            ,ab_mag_gaia=sp.vsep(uvstar_noab,uvstar_ab)
            ))
  raise
