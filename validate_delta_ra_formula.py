"""
validate_delta_ra_formula.py

Validate formula that calculates half-RA (Right Ascension) difference
of two planes that contain a conical FOV (Field Of View) with two known
parameters:
- Cone half-angle
- Declination (Dec) of the cone's axis.

Usage
=====
  [DEBUG=] python validate_delta_ra_formula.py [--hang=30.0]

Output
======
  (Declination,(half-RA calculate, half-RA estimate)) rows, to stdout

Derivation
==========
  See triple-quoted string [validate_delta_ra_formula.py.derivation] at
  the end of this script file

"""

import os
import sys
import math
import spiceypy as sp
import traceback as tb

do_debug = 'DEBUG' in os.environ

########################################################################
### Parse command-line argument(s), setup parameters

### --hang=<cone half-angle, degrees>
default_hang = 30.0
hangdeg = float(([default_hang] + [s[7:] for s in sys.argv[1:]
                                   if s.startswith('--hang=')
                                  ]
                         ).pop()
                        )

### Conversion factors between degrees and radians
dpr,rpd = sp.dpr(),sp.rpd()

### Cone half-angle and samples of Declinations, from 0 to (89.5-delta)
decdegs = [(0.5*decdec)-(decdec>0 and 1e-6 or 0) for decdec in range(180)]

### Convert to radians
hangrad = rpd * hangdeg
decrads = [rpd*decdeg for decdeg in decdegs]

### Trig functions of same
tansqhang = math.tan(hangrad)**2
sinhang = math.sin(hangrad)
coshang = math.cos(hangrad)

### List of triples of tan(Dec) squared, cos(Dec), and Dec
trigdecs = [(math.tan(decrad)**2,math.cos(decrad),decrad) for decrad in decrads]

### Unit basis vectors used in estimation of half-RA
vx,vy,vz = sp.ident()
### - Rotation X axis by cone half-angle around Z
vcone0 = sp.vrotv(vx,vz,hangrad)
### - Rotate vcone around X axis to represent circumference of cone's
###   intersection with celestial sphere at Radius=1; these vectors will
##    be rotated by angle -Declination around +Y axis, and their RA
##    calculated to estimate the half-RA at that declination
##    - Use cosine to oversample rotations near 0 and PI/2
vcones = [sp.vrotv(vcone0,vx,sp.halfpi()*(1+math.cos(sp.pi()*i/1e3))/2) for i in range(1000)]

########################################################################
### Initialize output list, loop over declinations
halfras = list()
for tansqdec,cosdec,decrad in trigdecs:
  try:

    ### Apply formula at this Dec(lination) and cone half-angle (hang)
    ### to get half-RA
    ###   arctan( sqrt(1-(tan(Dec)*tan(hang)**2)) / (cos(dec)*cos(hang))
    T = sinhang / math.sqrt(1.0 - (tansqdec * tansqhang))
    halfracalc = dpr * math.atan(T / (cosdec * coshang))

    ### Make estimate of same via the maximum RA of many cone vectors
    ### rotated by -Dec around Y
    halfraest = dpr * max([sp.recrad(sp.vrotv(vcone,vy,-decrad))[1]
                            for vcone in vcones
                           ])

    ### Append results to output list
    halfras.append((halfracalc-halfraest,halfracalc,halfraest,))

  except ValueError as e:
    ### At the point when the sum of (Dec + cone-half-angle) is 90 or
    ### more, the sqrt argument will be negative; confirm that is the
    ### case and exit the loop
    if do_debug: tb.print_exc()
    assert 1e-12 > (sp.halfpi() - (decrad+hangrad))
    break
  except:
    ### Re-raise any other exception
    raise

########################################################################
### Truncate input declinations list, output results
decs = decdegs[:len(halfras)]
for tup in zip(decs,halfras): print(tup)

########################################################################
### Plot results
import matplotlib.pyplot as plt

fig,(ax0,ax1,) = plt.subplots(2,1,sharex=True)

fig.suptitle('Cone half-angle = {0}deg'.format(hangdeg))

ax0.set_ylabel('$\Delta$RA, deg')
ax1.set_ylabel('Estimate error, microdeg')
ax1.set_xlabel('Declination, deg')

diffs,ayes,bees = zip(*halfras)
ax0.plot(decs,[b for b in bees],'o',label='Est.')
ax0.plot(decs,[a for a in ayes],label='Calc')
ax1.plot(decs,[diff*1e-6 for diff in diffs],'o',label='Calc-Est.')

ax0.legend(loc='upper left')
ax1.legend(loc='upper left')

plt.show()

derivation = """

Derivation
==========
  FOV cone intersecting a unit sphere is modeled as an axis
  perpendicular to, and ending at, a disk on the end, with the disk
  circumference on the unit sphere.  The cone half-angle is [hang]; the
  radius of the disk is sin(hang); the distance along the cone axis from
  the center of the sphere to the disk is cos(hang).

  We define [RA-constant planes] as planes that contain the sphere's
  polar axis.

  We define a point [T] to be the point along the line of the disk
  diameter, which diameter is parallel to the equator, where the space
  between an [RA-constant plane containing [T] and one point of the
  disk] and an [RA-constant plane containing the cone axis] contains
  exactly half of the disk.  The angle between those two planes is the
  angle that is to be solved for.

  Notation:

    [P]    -  Point "P"
    [PQ]   -  line from [P] to [Q]
    |PQ|   -  length pf [PQ]
    [PQR]  -  Angle from [P] to [Q] (vertex) to [R]; alsso triangle

  When the cone axis is on the equator, Dec(lination) = 0, and, looking
  down from the top (from +Z), the view projected onto the equatorial
  plane is like this:


                        -->|     |<--sin(hang) = T
                           |     |
       Unit sphere       | |     |
                 |       v    __...__
                 v   ___.--T=====@=====T--.__  <-disk viewed edge-on
                _.--'       \    |           `--._
             _-'             \_.-|<-half-RA = hang`-_
            /                 \  |                   \
           /                   \ |                    \
          |                     \|                     |
         |              -------  o <-center of sphere   |
                         ^       ^    also polar axis
                         |       |
                 cos(hang)       cone axis = [@o]


  In this case, [T] is at the end of the disk diameter parallel to the
  equator, and the half-RA angle, i.e. angle [@oT], is the same as hang,
  the half-angle of the cone.

  As the axis moves off the equator, i.e. to non-zero Dec values, the
  circumference of the disk will project onto the equatorial plane as an
  ellipse:

                        -->|     |<--sin(hang)
                           |     |
                           |     |
      sin(hang)sin(Dec)    |     |
                      |    |  -->|         |<--To-be-solved T
  cos(hang)cos(Dec)   |     ____...____    |
                  |  _|_.--'__..---..__`--. _
                _-v-' v    /<-ellipse->\     `--._
             _-' --------  |-----@-----| --T      `-_
            /              \__   |   __/_-'          \
           /         ----     ``---''_-'              \
          |           ^          |_-'                  |
         |       -----|--------  o                      |
                  ^
                  |


  In this case, the line [oT] is tangent to the projected ellipse, so
  when we solve for the horizontal distance to T, we can get the half-RA
  angle, i.e. angle [@oT], here as arctangent(T / (cos(hang)cos(Dec)),
  where the denominator is the projected length of the cone axis onto
  the equatorial plane.

  To solve for distance [@T], we stretch the projection along the
  projection of the cone axis by 1/sin(Dec), so the ellipse will project
  as a circle of radius sin(hang), and the line [oT] will be tangent to
  that circle at point [#] in the diagram below


                      sin(hang)->|     |<--
                                 |     |
                                 |     |
               sin(hang)         |     |
                       |      -->|     |   |<--To-be-solved [@T]
  cos(hang)/tan(Dec)   |         |     |   |
                   |   |         |     |   |
                   v   v
                  -------------  @-----| --T
                                 |\_   | _/
                                 |  \_ //
                                 |  _.#/
                     ----------  --' /
                      ^          | _/
                      |          |/
                  -------------  o
                   ^
                   |

  Since [oT] is tangent at [#] and the ellipse is now a circle, the
  angle [o#@] is a right angle, triangle [o#@] is a right triangle with
  [o@] as the hypotenuse, and the leg [@#] is a radius of the circle so
  its length is sin(hang), so length

            _______________________________________
    |o#| = v (cos(hang)/tan(Dec))**2 - sin(hang)**2

  Angle [o@T] is a right angle and a right triangle, as is [@#T], and
  all three triangles are similar to each other, with angles [@o#],
  [To@], and [T@#] being equal.  From that, we know the following ratios
  are equal:

    |@T|     |o@|                     |@#| |o@|
    ----  =  ----  =>  |@T|  =  T  =  ---------
    |@#|     |o#|                       |o#|


               sin(hang)  cos(hang) / tan(Dec)
    T  =  -----------------------------------------
           ________________________________________
          v (cos(hang)/tan(Dec))**2 - sin(hang)**2


                        sin(hang)
    T  =  -----------------------------------------
               ____________________________
              v 1 - (tan(hang) tan(Dec))**2


   And finally, because T in the stretched projection is the same as T
   from the original un-stretched projection:

                      ( cos(Dec) cos(hang) )
     halfRA   = arctan( ------------------ )
                      (         T          )

"""
