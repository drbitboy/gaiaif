"""
Transformation between J2000 to ICRS inertial reference frames [0]


Summary
=======

TK data expressing the following rotations [1]

  J2000 Pole vector, +Z, wrt ICRS Pole vector:

  1)  6.8mas toward 18h => + 6.8mas around +X => rotate frame - 6.8mas around +X
  2) 16.6mas toward 12h => -16.6mas around +Y => rotate frame +16.6mas around +Y

  J2000 Equinox wrt ICRS Equinox:

  3) 78.0mas toward 12h => +78.0mas around +Z => rotate frame -78.0mas around +Z


Usage
=====

- SPICE Kernel:

      DOUBLEPRECISION MTX_ICRS_TO_J2000(3,3)
      ...
      CALL FURNSH('gaia_icrs_tk.py')
      CALL PXFORM('ICRS','J2000',0.0,MTX_ICRS_TO_J2000)

- Python script (BASH command line, with % as prompt):

      % [DEBUG=] python gaia_icrs_fk.py && echo Success || echo Failed


Data [1]
========

TK Frame as a SPICE Frames Kernel (FK) [1]

\begindata
FRAME_ICRS               = 1700111
FRAME_1700111_NAME       = 'ICRS'
FRAME_1700111_CLASS      = 4
FRAME_1700111_CLASS_ID   = 1700111
FRAME_1700111_CENTER     = 10

TKFRAME_1700111_RELATIVE = 'J2000'
TKFRAME_1700111_SPEC     = 'ANGLES'
TKFRAME_1700111_UNITS    = 'ARCSECONDS'
TKFRAME_1700111_AXES     = ( 1, 2, 3 )
TKFRAME_1700111_ANGLES   = ( +6.8D-3
                           , -16.6D-3
                           , +78.0D-3
                           )
\begintext


References
==========

[0] SPICE FK Required Reading

[1] M. Gontier, E.F. Arias, and C. Barache, "Maintenance of the ICRF
    using the most stable sources," ICRS Center Report for 2001-2004,
    (IERS Technical Note No. 34), Jean Souchay and Martine
    Feissel-Vernier (eds.), p. 9, Fig. 1; see figure in ASCII below.
==> http://www.iers.org/MainDisp.csl?pid=46-1100077
==> https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn34.html
==> https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote34/tn34.pdf?__blob=publicationFile&v=1

N.B. there is also an earlier publication frm the US Naval Observatory
that has slightly different numbers:  5.1mas and 17.3mas around X and Y;
also 78mas around Z.  Cf.  McCarthy, D.D., "IERS Conventions (1992),"
IERS Technical Note 21, July 1996,
=> https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote21/tn21.pdf?__blob=publicationFile&v=1


Figure 1 from [1]; P is pole (+Z); origin is ICRS pole 

                                 RA=12h
                                   |
                   P               |
                    J2000          |
                          o-----------
                          |        | ^
                          |        | 16.6mas toward 12h
                          |        |
                          |        |
                          |        |
                          |        |
              |<19.9mas   |        |
              | toward    |        |
              | 18h       |        |
  RA=18h------|-----------|--------+------------------>RA=06h
              |           |<6.8mas |P                  +Y
              |             toward | ICRS
              |             18h    |
              |                    |
              |                    | 9.1mas toward 00h
              |                    | v
          P   o-----------------------
           FK5                     |
                                   V
                                 RA=00h
                                   +X

Equinoxes (+X)
- O = origin of Right Ascensions (RAs) of ICRS and FK5
- A = First point of Ares i.e. Spring equinos of J2000

          |             |                      |
          |<--22.9mas-->|<-------78mas-------->|
          |             |                      |
  --------o-------------+----------------------o-------------> +RA
         O             O                      A
          FK5           ICRS                   J2000

                           


Python script
=============

The rest of the lines in this file, after this docsstring, conmprise a
test of these data, written in Python; https://www.python.org/

"""
import os
import sys
import spiceypy as sp
from astropy.coordinates import SkyCoord

do_debug = 'DEBUG' in os.environ

### Conversion factor:  milliarcsecond/radian = mas/deg * deg/rad
maspr = 3.6e6 * sp.dpr()


########################################################################
### I.  Test the methodology
########################################################################

### Astropy provides FK5 Pole and Equinox vectors in ICRS frame;
### SPICE TWOVEC creates rotation matrix from those two vectors
fk5_z_in_icrs_astropy = SkyCoord(frame='fk5',x=0,y=0,z=1,representation_type='cartesian').icrs.cartesian.xyz.value
fk5_x_in_icrs_astropy = SkyCoord(frame='fk5',x=1,y=0,z=0,representation_type='cartesian').icrs.cartesian.xyz.value
icrs_to_fk5_astropy = sp.twovec(fk5_z_in_icrs_astropy,3,fk5_x_in_icrs_astropy,1)
###
### FK5 Pole vector, +Z, wrt ICRS Pole vector:
### 1) 19.9mas toward 18h => +19.9mas around +X => rotate frame -19.9mas around +X
### 2)  9.1mas toward 00h => + 9.1mas around +Y => rotate frame - 9.1mas around +Y
### FK5 Equinox wrt ICRS Equinox:
### 3) 22.9mas toward 180 => -22.9mas around +Z => rotate frame +22.9mas around +Z
### SPICE MXM and XPOSE creates rotation matrix combining those rotations
fk5_to_icrs_199x = sp.rotate(-19.9/maspr,1)
fk5_to_icrs_091y = sp.rotate(- 9.1/maspr,2)
fk5_to_icrs_229z = sp.rotate(+22.9/maspr,3)
icrs_to_fk5_rots = sp.xpose(sp.mxm(sp.mxm(fk5_to_icrs_229z,fk5_to_icrs_091y),fk5_to_icrs_199x))

assert 2e-14>abs((icrs_to_fk5_rots-icrs_to_fk5_astropy)/icrs_to_fk5_astropy).max()


########################################################################
### II. Use same methodology to generate ICRS to J2000 rotation matrix
########################################################################

j2k_to_icrs_068x = sp.rotate(- 6.8/maspr,1)
j2k_to_icrs_166y = sp.rotate(+16.6/maspr,2)
j2k_to_icrs_780z = sp.rotate(-78.0/maspr,3)
icrs_to_j2k_rots = sp.xpose(sp.mxm(sp.mxm(j2k_to_icrs_780z,j2k_to_icrs_166y),j2k_to_icrs_068x))


########################################################################
### III. Get the transform from the FK data in the docstring above
########################################################################

sp.furnsh(__file__)
icrs_to_j2k_fk = sp.pxform('icrs','j2000',0.0)
sp.unload(__file__)

### Calculate the difference
frac_diff=(icrs_to_j2k_fk-icrs_to_j2k_rots) / icrs_to_j2k_rots

### Output if requested
if do_debug:
  import pprint
  pprint.pprint(dict(icrs_to_j2k_rots=icrs_to_j2k_rots*maspr
            ,icrs_to_j2k_fk=icrs_to_j2k_fk*maspr
            ,j2k_x_in_icrs=sp.mtxv(icrs_to_j2k_rots,[1,0,0])*maspr
            ,j2k_z_in_icrs=sp.mtxv(icrs_to_j2k_rots,[0,0,1])*maspr
            ,frac_diff=frac_diff
            )
     )

### Test the maximum difference
assert 2e-14>abs(frac_diff).max()
