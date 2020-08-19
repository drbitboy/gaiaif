import os
import sys
import math
import numpy
import spiceypy as sp

try: dpr
except:
  dpr = sp.dpr()
  rpd = sp.rpd()
  rpmas = rpd / 3600e3              ### Radian / milliarcsecond
  aupkm = sp.convrt(1.,'km','au')   ### Astonomical Unit / kilometer
  recip_clight = 1.0 / sp.clight()


########################################################################
class FOV(object):
  """
Field-Of-View (FOV) and methods to determine if a ray is inside the FOV

Attributes

.L            Length of FOV argument to .__init__
.fovtype      FOV type:  FOV.CIRCLETYPE; RADECBOXTYPE;,.POLYGONTYPE
.radecdegs    List of [RA,Dec] pairs of input vectors/vertices, degrees
.uvfovxyzs    List of Unit vectors of [X,Y,Z] triples of input vectors
.hangdeg      Cone half-angle for circle FOV, degrees
.hangrad      Cone half-angle for circle FOV, radians
.radec_boxes  List of lists:  FOV bounding boxes; ralo,rahi,declo,dechi
.convex       True is FOV is convex, else False

- FOV.POLYGONTYPE attributes:

.mtxtofov   Rotation matrix from inertial frame to local reference frame
.uvlclxyzs  Vertices in local reference frame (reffrm), unit vectors
.uvavg      Mean vector of all vertices, will be +Z of local reffrm
.localxyzs  Verts in local reffrm, on plane Z=+1, w/stellar aberration

"""
  (CIRCLETYPE,RADECBOXTYPE,POLYGONTYPE
  ,) = 'circle radecbox polygon'.split()

  ######################################################################
  def __init__(self,fovraws,obs_pos=None,obs_vel=None,obs_year=None):
    """Convert FOV definition in external inertial reference frame to an
FOV definition in a local reference frame; also determine RA,Dec limits
of FOV for use in star catalog lookup.

For polygonal FOVs, set up FOV as vectors in a local reference frame,
with a matrix to rotate vectors from the external to the local frame.
The local reference frame (reffrm) will have +Z along the average
vertices' direction, and will be rotated such the the +X axis is not
parallel to any of the edges of the polygon from the FOV projected onto
the Z=+1 plane.

Arguments

  fovraws - sequence either of vector and cone and half-angle for
            circular FOV, or of a pair of vectors for RA,Dec box, or
            of three or more vectors for a polygonal FOV.  A vector is
            a sequence either of two values, RA,Dec, or of three values,
            X,Y,Z.

  obs_pos - Observer position, solar system barycentric, 3-vector, km
            - For parallax correction

  obs_vel - Observer velocity, solar system barycentric, 3-vector, km/s
            - For proper motion correction

  obs_year - Observer time, y past 2015.5 (Gaia DR2 epoch)
            - For stellar aberration correction

"""
    ### Get count of items in argument sequence; ensure it is 2 or more
    (self.fovraws
    ,self.obs_pos
    ,self.obs_vel
    ,self.obs_year
    ,)= fovraws,obs_pos,obs_vel,obs_year
    self.L = len(fovraws)
    assert 1<self.L, 'Too few vertices in FOV'

    ################################
    ### Initialize:  FOV RA,Dec pairs; FOV type (assume polygon); FOV
    ###              vector triples; list of RA,Dec boxes
    self.radecdegs,self.fovtype = list(),FOV.POLYGONTYPE
    self.uvfovxyzs,fovsum = list(),sp.vpack(0.,0.,0.)
    self.radec_boxes = list()
    rdba = self.radec_boxes.append   ### Shorthand to append box to list

    ################################
    ### Parse list of vertices:
    ### - [list,float]          =>  Circle (cone)
    ### - [list,list]           =>  RA,Dec box
    ### - [list,list,list,...]  =>  Polygon
    for vertex in fovraws:

      ### For second of two vertices ...
      if 1==len(self.radecdegs) and 2==self.L:
        ### Two-vertex items are either a conic FOV, or an [RA,Dec] box
        try:
          ### If second item in list is a float, then it's a half-angle
          ### of the cone
          self.hangdeg = float(vertex)
          assert self.hangdeg < 90.0,'Cone half-angle is not less than 90degrees'
          assert self.hangdeg > 0.0,'Cone half-angle is not greater than 0degrees'
          self.hangrad = self.hangdeg * rpd
          self.min_cosine = math.cos(self.hangrad)
          self.uv_cone_axis = self.uvfovxyzs[0]
          self.fovtype = FOV.CIRCLETYPE
          break
        except AssertionError as e:
          raise
        except:
          ### If the above fails, then it's the second corner of the box
          self.fovtype = FOV.RADECBOXTYPE

      ### Parse one vertex
      ra,dec,uvxyz = parse_inertial(vertex)

      ### Append RA,Dec and unit vector XYZ onto their resepective lists
      self.radecdegs.append((ra,dec,))
      self.uvfovxyzs.append(uvxyz)
      fovsum = sp.vadd(fovsum,uvxyz)

    ################################
    ### Calculate RA,DEC limits as list of [ralo,rahi,declo,dechi] boxes
    ### - .radec_boxes is a list; rdba is .radec_boxes.append
    ### - List will have multiple RA,Dec boxes if FOV crosses the Prime
    ###   Meridian (PM) an even number of times.

    if self.fovtype == FOV.RADECBOXTYPE:
      ### RA,DEC box FOV:  calculate limits; handle PM crossing
      ras,decs = zip(*self.radecdegs)
      ralo,rahi = sorted(ras)
      declo,dechi = sorted(decs)
      if 180 > (rahi-ralo):
        rdba([ralo,rahi,declo,dechi])
      else:
        rdba([0.0,ralo,declo,dechi])
        rdba([rahi,360.0,declo,dechi])

    elif self.fovtype == FOV.CIRCLETYPE:
      ### Circular FOV:  DEC limits determine RA limits; handle PM Xing
      ra,dec = self.radecdegs[0]
      fovdeclo = dec - self.hangdeg
      fovdechi = dec + self.hangdeg

      if fovdeclo < -90.0 or fovdechi > 90.0:
        ### A pole is in the FOV; use full RA range
        fovralo,fovrahi = 0.0,360.0
        fovdeclo,fovdechi = max([fovdeclo,-90.0]),min([fovdechi,+90.0])

      elif fovdeclo == -90.0 or fovdechi == 90.0:
        ### A pole is on the FOV circumference; RA range is 180 degrees
        fovralo,fovrahi = ra-90.0,ra+90.0

      else:
        ### The FOV excludes the poles; calculate the RA range, using
        ### the formula validated in script validate_delta_ra_formula.py
        tanhang,tandec = math.tan(self.hangrad),math.tan(dec*rpd)
        sinhang,cosdec = math.sin(self.hangrad),math.cos(dec*rpd)
        coshang = math.cos(self.hangrad)
        T = sinhang / math.sqrt(1.0 - ((tanhang*tandec)**2))
        deltara = dpr * math.atan(T / (cosdec * coshang))
        fovralo,fovrahi = ra-deltara,ra+deltara

      ### Ensure RA limits are within range [0:360] (N.B. inclusive)
      if fovralo < 0.0: fovralo += 360.0
      if fovrahi > 360.0: fovrahi -= 360.0

      if fovralo <= fovrahi:
        ### RA lo <= RA hi:  no PM crosssing
        rdba([fovralo,fovrahi,fovdeclo,fovdechi])
      else:
        ### RA hi < RA hi:  there is a PM crosssing
        rdba([0.0,fovrahi,fovdeclo,fovdechi])
        rdba([fovralo,360.,fovdeclo,fovdechi])

    else:
      assert self.fovtype == FOV.POLYGONTYPE
      ### Polygonal FOV:  build frame where all vertices will be
      ### projected onto the plane Z=1

      ### .uvavg:  unit vector = mean of all vertices, will be +Z
      self.uvavg = sp.vhat(fovsum)

      ### Create rotation matrix to FOV frame:  +Z is mean of vertices'
      ###   directions (.uvavg); +X will be a direction that is not
      ###   parallel to any side of the polygon
      ### - Start with temporary matrix with +Z as defined above; +X
      ###   toward vertex at largest angle from .uvavg
      vother = min([(sp.vdot(self.uvavg,v),list(v),) for v in self.uvfovxyzs])[1]
      tmpmtx = sp.twovec(self.uvavg,3,vother,1)
      ### - Rotate all vectors to that frame; scale Z components to 1.0
      vtmps = list()
      for v in self.uvfovxyzs:
        ### - Ensure all vertices are in the same hemisphere
        assert 0.0 < sp.vdot(self.uvavg,v),'All vertices are not in the same hemisphere'
        vtmp = sp.mxv(tmpmtx,v)
        vtmps.append(sp.vscl(1.0/vtmp[2],vtmp))

      ### Find largest azimuth gap between any two sides:  that azimuth
      ###   will be direction of +X in the final rotation matrix
      ### - Get azimuths of all sides of polygon, in range [-PI:PI]
      azimuths,vlast = list(),vtmps[-1]
      for v in self.uvfovxyzs:
        azimuths.append(numpy.arctan((v[1]-vlast[1])/(v[0]-vlast[0])))
        vlast = v
      ### - Sort angles and add [least angle plus PI] to end of list
      azimuths.sort()
      azimuths.append(azimuths[0]+sp.pi())
      ### - Find largest delta-azimuth and its index
      dazimuths = [hi-lo for hi,lo in zip(azimuths[1:],azimuths[:-1])]
      maxdaz = max(dazimuths)
      imaxdaz = dazimuths.index(maxdaz)
      ### - Calculate azimuth from to mean of that delta-azimuth,
      meanaz = azimuths[imaxdaz] + (maxdaz / 2.0)

      ### Final matrix:  add rotation of tmpmtx around +Z by that angle
      self.mtxtofov = sp.mxm(sp.rotate(meanaz,3),tmpmtx)

      ### Apply final rotation matrix, store results in .uvlclxyzs
      tmpmtx = sp.twovec(self.uvavg,3,vother,1)
      self.uvlclxyzs = [self.rotate_to_local(v) for v in self.uvfovxyzs]

      ### Calculate upper and lower RA and Dec limits, with PM crossings
      los,his = list(),list()
      ### - Create [[RA,Dec],[X,Y,Z]] pairs list; ensure last is off PM
      pairs = list(zip(self.radecdegs,self.uvfovxyzs))
      pop_count = 0
      while pairs[-1][0][0] == 0.0:
        pop_count += 1
        assert pop_count < self.L,'All vertices are on the Prime Meridian'
        pairs.append(pairs.pop(0))

      ### Count PM crossings
      self.crossing_count = 0
      lastra = pairs[-1][0][0]
      zero_count = 0
      for (ra,dec,),xyz in pairs:
        if ra == 0.0:
          zero_count += 1
          if lastra > 180.0: ra = 360.0
        if 180 < abs(ra-lastra): self.crossing_count += 1
        lastra = ra

      if 0==self.crossing_count or 1==(1&self.crossing_count):
        ### If there are either no, or an odd number, of PM crossings,
        ### then use the pairs as-is for a single FOV
        subfovs = [pairs]
        if self.crossing_count:
          ### - For odd crossing count, one pole or the other must be
          ###   in the FOV; init full RA range, that pole for Dec ranges
          ralo,rahi = 0.0,360.0
          if sp.vdot(self.uvavg,[0,0,1]) > 0.0: declo = dechi = +90.0
          else                                : declo = dechi = -90.0
        else:
          ### - For zero crossing count, initialize inverted ranges
          ralo,rahi = 360.0,0.0
          declo,dechi = +90.0,-90.0
        subranges = [[ralo,rahi,declo,dechi]]

      else:
        ### If there are an even, non-zero number of PM crossings, break
        ### them into two sub-FOVs, one on either side of the PM

        eastfov,westfov = list(),list()

        if zero_count:
          ### If there are any zero RA values, rotate the pairs to
          ### ensure a zero-RA pair is the first, so it and the non-zero
          ### last pair will be assigned to the correct side of the PM
          while pairs[0][0][0]!=0.0: pairs.append(pairs.pop(0))
        else:
          ### If there are no zero RA values, rotate the pairs to ensure
          ### a crossing occurs between the last and first pair, so the
          ### corresponding zero crossing will be assigned to the
          ### correct side of the PM
          while abs(pairs[0][0][0]-pairs[-1][0][0])<180:
            pairs.append(pairs.pop(0))

        ### Write vertices into the two sub-FOVs

        ### - Set last-vertex values for first item in pairs
        (lastra,lastdec,),lastxyz = pairs[-1]

        for pair in pairs:
          ### - Loop over vertex pairs ((RA,DEC,),Cartesian_Vector)
          (ra,dec,),xyz = pair

          if ra == 0.0:

            ### - When RA=0, the previous RA determines if it's 0 ar 360
            if lastra >= 180.0:
              ra = 360.0
              westfov.append([(ra,dec,),xyz])
              iswest = True
            else:
              eastfov.append(pair)
              iswest = False

          elif abs(lastra-ra) >= 180.0:

            ### - When the change in RA>=180, the PM is being crossed

            ### - Find the mid-vector where the PM is crossed
            k1 = -xyz[1] / (lastxyz[1]-xyz[1])
            midxyz = sp.vhat(sp.vlcom(1.0-k1,xyz,k1,lastxyz))
            middec = dpr * sp.recrad(midxyz)[2]

            ### - Add that mid-vector, with RA=360, to the west FOV
            westfov.append([(360.0,middec,),midxyz])

            ### - Determine if vector is west
            iswest = ra >= 180.0

            ### - Add that mid-vector, with RA=0, to the east FOV ...
            if (ra > 0.0) and (not iswest):
              ### - ... only if the ra is not already 0, as it will be
              ###       added in the next step
              eastfov.append([(0.0,middec,),midxyz])

            ### Add the vector to either east or west FOV
            if iswest: westfov.append(pair)
            else     : eastfov.append(pair)

          else:

            ### PM was not crossed, add vector to same FOV, as last time
            if iswest: westfov.append(pair)
            else     : eastfov.append(pair)

          ### - Set last-vertex values for next item in pairs
          (lastra,lastdec,),lastxyz = (ra,dec,),xyz

        ### - Create subfovs list of east and west FOVs; set subranges
        subfovs = [eastfov,westfov]
        subranges = [[360.0,0.0,90.0,-90.0],[360.0,0.0,90.0,-90.0]]

      ### To here, we have list of FOV(s) and list of range(s); use them
      ### to determine RA,DEC box(es) to use for database query

      while subfovs:

        ### Get sub-FOV, sub-range; set last vertex's XYZ
        subfov,(ralo,rahi,declo,dechi,) = subfovs.pop(),subranges.pop()
        lastxyz = subfov[-1][-1]

        for pair in subfov:
          ### Each element of subfov comprises (RA,Dec) and vertex XYZ
          ### - xyz is a unit vector
          (ra,dec,),xyz = pair

          ### - Adjust RA limits as needed from RA of vertex
          if   ra > rahi: rahi = ra
          elif ra < ralo: ralo = ra

          ### - Set Dec extrema from DEC of vertex
          maxdec = mindec = dec

          ### - Calculate Dec extrema from lastxyz to xyz
          ### -- Normal to plane of lastxyz and syz
          sidenormal = sp.vcrss(lastxyz,xyz)
          ### -- Z-rates along great circle at lastxyz and at xyz
          lastdz = sp.vcrss(sidenormal,lastxyz)[2]
          dz = sp.vcrss(sidenormal,xyz)[2]
          if 0.0 > (lastdz*dz):
            ### -- If sign of Z-rates differs, there should be an
            ###    extreme value between lastxyz and xyz
            ### --- Get vector perpendicular to side normal on equator
            ### --- Use that to calculate the unit vector at Dec extreme
            equinox = sp.vcrss([0,0,1],sidenormal)
            vtoextremez = sp.ucrss(sidenormal,equinox)
            ### --- Cosine of angle between lastxyz and xyz
            mindot = sp.vdot(lastxyz,xyz)
            for none in [None,None]:
              ### --- Two cases:  vtoextremez and -vtoextremez
              ###     - Angles from vtoextremez to lastxyz and to xyz
              ###       must be less than angle between lastxyz and xyz
              ###       so cosines of those angles must be greater
              lastxyzdot = sp.vdot(lastxyz,vtoextremez)
              xyzdot = sp.vdot(xyz,vtoextremez)
              if lastxyzdot>mindot and xyzdot>mindot:
                ### --- Adjust maxdec and mindec as needed
                try   : extremedec = dpr * math.asin(vtoextremez[2])
                except: extremedec = dpr * sp.recrad(vtoextremez)[2]
                if   extremedec > maxdec: maxdec = extremedec
                elif extremedec < mindec: mindec = extremedec
                break
              ### --- Invert vtoextremez for next pass
              vtoextremez = sp.vminus(vtoextremez)

          ### - Adjust Dec limits as needed from Dec extrema of side
          if maxdec > dechi: dechi = maxdec
          if mindec < declo: declo = mindec
          lastxyz = xyz

        ### Append calculated RA,Dec box(es)
        rdba((ralo,rahi,declo,dechi,))

      ### Put None in .localxyzs, in .v_for_stellar_aberr, and in
      ### .v_for_parallax; if no stellar aberration or parallax is
      ### explicitly applied to define it later, then .localxyzs will be
      ### calculated on the fly
      self.localxyzs = None
      self.v_for_stellar_aberr = None
      self.v_for_parallax = None

  ########################################################################
  def __repr__(self): return str(self.fovraws)

  ########################################################################
  def star_in_fov(self
                 ,vstar
                 ,parallax_maspau=None
                 ,pmra_maspy=None
                 ,pmdec_maspy=None
                 ):
    """Return True if the star (RA,Dec or xyz) argument is in the FOV

Argument vstar is either an RA,Dec pair (degrees) or a 3-vector

"""

    if self.fovtype == FOV.RADECBOXTYPE:
      ##################################################################
      ### Compare star vector to [RA,Dec] box
      ra,dec = parse_inertial(vstar,return_radec=True)[:2]
      for ralo,rahi,declo,dechi in self.radec_boxes:
        if ra<ralo: continue
        if ra>rahi: continue
        if dec<declo: continue
        if dec>dechi: continue
        return True,sp.radrec(1.,rpd*ra,rpd*dec)
      return False,None

    ### Get inertial star unit vector without RA,Dec
    uvinertial = uvraw = parse_inertial(vstar,return_radec=False)

    ### Corrections for direction to star
    ### - Assume all corrections are small and can be applied in units
    ###   of radians to a unit vector

    ### - Proper Motion
    if (not (None is self.obs_year)
       ) and (not (None in (pmra_maspy,pmdec_maspy,))
       ) and (pmra_maspy != 0.0 or pmdec_maspy != 0.0):
      uveast = sp.ucrss([0,0,1],uvraw)
      uvnorth = sp.ucrss(uvraw,uveast)
      ###cosdec = math.sqrt(1.0 - (uvraw[2]*uvraw[2]))
      uvinertial = sp.vhat(sp.vlcom3(self.obs_year*rpmas*pmdec_maspy,uvnorth
                                    ,self.obs_year*rpmas*pmra_maspy,uveast
                                    ,1.0,uvinertial
                                    )
                          )

    ### - Parallax
    if (not (None is self.obs_pos)
       ) and (not (None is parallax_maspau)
       ) and (parallax_maspau != 0.0):
      uvinertial = sp.vhat(sp.vsub(uvinertial
                                  ,sp.vscl(aupkm*parallax_maspau*rpmas,self.obs_pos)
                                  )
                          )

    ### - Stellar Aberration
    if not (None is self.obs_vel):
      uvinertial = sp.vhat(sp.vadd(uvinertial
                                  ,sp.vscl(recip_clight,self.obs_vel)
                                  )
                          )

    if self.fovtype == FOV.CIRCLETYPE:
      ##################################################################
      ### Compare inertial star vector to circular FOV
      return sp.vdot(uvinertial,self.uv_cone_axis) >= self.min_cosine,uvinertial

    assert FOV.POLYGONTYPE == self.fovtype,'Unknown FOV type [{0}]'.format(self.fovtype)

    ####################################################################
    ### Compare star vector to polygonal FOV

    if self.is_convex():
      ### Convex FOV:  a negative dot product with the inward-pointing
      ###              normal to any side indicates star is outside FOV
      for inwardnorm in self.inwardsidenorms:
        if sp.vdot(uvinertial,inwardnorm) < 0.0: return False,uvinertial

      ### All dot products were non-negative:  star is within FOV
      return True,uvinertial

    ### Rotate inertial unit vector to local reference frame (reffrm)
    uvlocalstar = self.rotate_to_local(uvinertial)

    ### Scale to Z=unity
    if uvlocalstar[2] < 1e-15: return False,uvinertial
    z1star = sp.vscl(1.0/uvlocalstar[2],uvlocalstar)

    ### Setup .localxyzs and .fovsides
    if None is self.localxyzs: self.setup_localxyzs()

    ### Count number of crossings of FOV sides
    count = len([None
                 for fovside in self.fovsides
                 if fovside.right_of(z1star)
                ])

    return ((count&1) and True or False),uvinertial


  ########################################################################
  def is_convex(self):

    try:
      assert hasattr(self,'convex')
      return self.convex
    except: pass

    assert self.fovtype == FOV.POLYGONTYPE, 'FOV.is_convex() method must called from polygon, not from {0}, FOV'.format(self.fovtype)

    uvuseds = list()
    uvlefts = self.uvfovxyzs[::]
    oneneg,onepos = False,False
    self.inwardsidenorms = list()

    norm = sp.ucrss(uvlefts[-1],uvlefts[0])
    for uv in uvlefts[1:-1]:
      dot = sp.vdot(norm,uv)
      if   0.0<dot: onepos = True
      elif 0.0>dot: oneneg = True
      else:
        self.convex = False
        return False

    if onepos: self.inwardsidenorms.append(norm)
    else     : self.inwardsidenorms.append(sp.vminus(norm))

    uv0 = uvlefts.pop()
    while uvlefts and (not (oneneg and onepos)):

      uv1 = uv0
      uv0 = uvlefts.pop()

      norm = sp.ucrss(uv0,uv1)

      for uv in (uvuseds+uvlefts):
        dot = sp.vdot(norm,uv)
        if   0.0<dot: onepos = True
        elif 0.0>dot: oneneg = True
        else:
          self.convex = False
          return False

      if onepos: self.inwardsidenorms.append(norm)
      else     : self.inwardsidenorms.append(sp.vminus(norm))

    self.convex = onepos ^ oneneg
    return self.convex


  ########################################################################
  def setup_localxyzs(self):
    """Scale unit vectors to Z=unity; build FOV sides"""
    self.localxyzs = [sp.vscl(1.0/v[2],v) for v in self.uvlclxyzs]
    self.fovsides = list()
    lastxyz = self.localxyzs[-1]
    for xyz in self.localxyzs:
      self.fovsides.append(FOVSIDE(xyz,lastxyz))
      lastxyz = xyz

  ########################################################################
  def rotate_to_local(self,vxyz):
    """Utility to rotation from inertial frame to local frame"""
    return sp.mxv(self.mtxtofov,vxyz)

  ########################################################################
  def setup_stellar_aberration(self,observer_velocity_xyz):
    """Scale observer velocity in the inertial frame (J2000 or ICRS)
by the reciprocal of speed of light.  The resulting vector can be added
to unit vectors, or subtracted, to correct for stellar aberration.

N.B. Stellar aberration correction is not yet implemented

"""
    self.v_for_stellar_aberr = sp.vscl(recip_clight,observer_velocity_xyz)

  ########################################################################
  def setup_parallax(self,observer_position_xyz):
    """Store observer position in the intertial frame (J2000 or ICRS);
it will be used in parallax correction

N.B. Parallax correction is not yet implemented

"""
    self.v_for_parallax = sp.vpack(*observer_position_xyz)

  ########################################################################
  def get_radec_boxes(self):
    """Retrieve RA,Dex boxes e.g. for SQL queries of a star catalog
Return:  [[Box0_RAlo,Box0_RAhi,Box0_Declo,Box0_Dechi]]
    OR   [[Box0_RAlo,Box0_RAhi,Box0_Declo,Box0_Dechi],[Box1_RAhi,...]]

"""

    return self.radec_boxes[:]


########################################################################
class FOVSIDE(object):
  """Encapsulate one side of an FOV"""
  def __init__(self,v0,v1):
    """Two endpoints of an FOV side; it is assumed that the Z components
are unity.

Setup .m and .b for   xside = .m*y + .b

    """
    self.vinputs = v0,v1
    self.xhi = max([v0[0],v1[0]])
    self.yhi,self.ylo = v0[1]>v1[1] and (v0[1],v1[1],) or (v1[1],v0[1])

    self.m = (v1[0]-v0[0]) / (v1[1]-v0[1])  ### (x1-x0)/(y1-y0)
    self.b = v0[0] - (v0[1] * self.m)       ### x0 - y0*(x1-x0)/(y1-y0)

  def right_of(self,v):
    """Is this side to the right of the vector argument v, with Z=unity,
or does this side pass over the vector; exclude the maximum y of this
side.

"""
    x,y = v[0:2]
    if y <  self.ylo: return False
    if y >= self.yhi: return False
    if x >  self.xhi: return False
    if x > ((y * self.m) + self.b): return False
    return True

########################################################################
def parse_inertial(input_vertex, return_radec=True):
  """Parse input vertex as either RA,Dec or XYZ;
Return RA,DEC (optional) and unit vector

"""
  if 2==len(input_vertex):
    ### Vertex has two items:  assume they are RA and Dec
    ra,dec = map(float,input_vertex)
    assert ra<=360.0 and ra>=0.0,'RA ({0}) is out of range [0,360)'.format(ra)
    assert dec<=90.0 and dec>=-90.0,'RA ({0}) is out of range [-90,+90]'.format(dec)
    uvxyz = sp.radrec(1.0,ra*rpd,dec*rpd)
  else:
    ### Vertex has three items:  assume they are XYZ
    assert 3==len(input_vertex),'XYZ input vector [{0}] for vertex does not have 3 elements'.format(str(input_vertex))
    uvxyz = sp.vhat(sp.vpack(*map(float,input_vertex)))
    if return_radec: ra,dec = sp.vsclg(dpr,sp.recrad(uvxyz)[1:],2)

  if return_radec: return ra,dec,uvxyz
  return uvxyz
