import os
import sys
import math
import gaiaif_util
import spiceypy as sp

rpd = sp.rpd()

if "__main__" == __name__:

  ### RA,Dec box
  ### Circle
  ### Polygon (non-convex quadrilateral)
  ### Polygon (convex quadrilateral)
  for fov in (gaiaif_util.FOV([[315,89.2],[15,89.8]])
             ,gaiaif_util.FOV([[345,89.3],0.5])
             ,gaiaif_util.FOV([[30,89],[0,88.5],[5.0,88.9],[270,89]])
             ,gaiaif_util.FOV([[30,89],[0,88.5],[315,88.5],[270,89]])
             ,gaiaif_util.FOV([[225,89],[0,88.5],[315,88.5],[90,89]])
             ,):

    radec_boxes = fov.get_radec_boxes()
    ###for radec_box in radec_boxes: print(radec_box) 
    if fov.POLYGONTYPE == fov.fovtype:
      print(dict(fov_is_convex=fov.is_convex()))

    def inbox(ra,dec):
      for ralo,rahi,declo,dechi in radec_boxes:
        if ra<ralo: continue
        if ra>rahi: continue
        if dec<declo: continue
        if dec>dechi: continue
        return True
      return False

    count = 650
    countplus = count + 1

    dtheta = 360.0 / count
    dphi = 2.0 / count
    ras = [i*dtheta for i in range(count)]*countplus
    decs = sum([[88+(i*dphi)]*count for i in range(countplus)],[])
    rtn = [(fov.star_in_fov(radec),radec,) for radec in zip(ras,decs) if inbox(*radec)]

    print(len(rtn))

    if not ('--no-plot' in sys.argv[1:]):
      import matplotlib.pyplot as plt

      xs,ys,zs = zip(*[sp.radrec(1.,*sp.vsclg(rpd,radec,2)) for inside,radec in rtn if inside])

      plt.axhline(0,color='lightgray')
      plt.axvline(0,color='lightgray')

      if len(xs) < len(rtn):
        outxs,outys,outzs = zip(*[sp.radrec(1.,*sp.vsclg(rpd,radec,2)) for inside,radec in rtn if not inside])
        plt.plot([x/rpd for x in outxs],[y/rpd for y in outys],',r',label='Outside FOV')

      plt.plot([x/rpd for x in xs],[y/rpd for y in ys],'.g',markersize=.3,label='Inside FOV')
      plt.xlabel('Colatitude wrt PM, ~deg')
      plt.ylabel('Colatitude wrt RA=+90deg Meridian, ~deg')
      plt.legend()
      plt.title('{1} FOV = {0}'.format(fov,fov.fovtype.upper()))
    
      plt.show()

