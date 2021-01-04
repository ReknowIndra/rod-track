# -*- coding: utf-8 -*-
"""
Example of module rodtrack
expects 119 images in the format "rod0000x.png" in subdirectory rod

Created on Mon Jan  4 13:08:08 2021

@author: Erwin K. Reichel

"""

import matplotlib.pyplot as plt
from rodtrack import *

# load images
ir = rod(r"rod/rod{0:05d}.png",nmax=119)

# define tracking region
ir.rodrect(172,480,30,150,23.0*pi/180)
ir.rsshow(0)
plt.plot(ir.cc,ir.cr,'r+')

# tracking with consecutive updates
#ir.seqpath()

# tracking with original region
ir.seqapath()

# show tracking path
pt = ir.showpath()
#imsave("trackpath.png",pt)

#%% print results and export image sequenc

print ("xp = {}".format(ir.xp))
print ("yp = {}".format(ir.yp))
print ("rot = {}".format(ir.rs))

# export all images showing path
for k in range(len(ir.ims)):
    imshow(ir.ims[k])
    plt.plot(ir.xp[:k],ir.yp[:k],'r--',linewidth=3.0)
#    ptk = ir.showpath(k)
#    imsave("pt{0}.png".format(k),ptk)

#%% show tracking path in red
ptk = ir.showpath(-1,rectonly=True)
imshow(ptk)
plt.plot(ir.xp,ir.yp,'r--',linewidth=3.0)