# -*- coding: utf-8 -*-
"""
Tracking of rodlike particles in microscopic video sequences

Created on Mon Jan  4 12:56:23 2021
Updated on Thu Mar  4 11:37:51 2021

@author: Erwin K. Reichel, JKU Linz, Austria
"""


from numpy import * 
from skimage.io import *
import matplotlib.pyplot as plt
from scipy import signal
import scipy.optimize as opt
from skimage.draw import line_aa, line, polygon, polygon_perimeter
import random

class rod:
    """stores image sequence and performs rodtracking"""
    
    def __init__(self,fn,nmax=10,dep=8):
        """loads nmax images with filenames fn.format(k+1) with color depth dep"""
        self.ims = []
        for k in range(nmax):
            try:
                print ("reading: "+fn.format(k+1))
                self.ims.append(array(imread(fn.format(k+1)),dtype=float)/2**dep)
            except:
                if k == 0:
                    raise IOError("No image file found.")
                else:
                    print("Images read: {0:d}".format(k))
                    break

    def rodregion(self,rr,rc):
        """define region to be tracked"""
        r,c = polygon(rr,rc)
        self.rr = array(r,dtype=float)
        self.rc = array(c,dtype=float)
        self.cent()
    
    def cent(self):
        """calculates center of tracking regio"""
        self.cr = mean(self.rr)
        self.cc = mean(self.rc)
        self.dcr = self.rr-self.cr
        self.dcc = self.rc-self.cc

    def rodrect(self,x,y,w,l,al=0):
        """coordinates for rotated (angle al) rectangular region with center of rotation at x,y"""
        rr = array([y,y-w*sin(al),y-w*sin(al)+l*cos(al),y+l*cos(al)])
        cc = array([x,x+w*cos(al),x+w*cos(al)+l*sin(al),x+l*sin(al)])
        self.ser = rr
        self.sec = cc
        self.rodregion(rr,cc)
        self.x = x
        self.y = y
        self.w = w
        self.l = l
        self.al = al
        
    def rsshow(self,n):
        """shows image with polygon region in white (color code 1)"""
        im = copy(self.ims[n])
        spr,spc = polygon_perimeter(self.ser,self.sec)
        im[spr,spc] = 1
        imshow(im)

    def rodtransform(self,rot,trans):
        """rotate and translate the tracking region"""
        corot = cos(rot)
        sirot = sin(rot)

        dcc = self.dcc*corot-self.dcr*sirot+trans[0]
        dcr = self.dcc*sirot+self.dcr*corot+trans[1]

        self.rr = self.cr + dcr
        self.rc = self.cc + dcc
        self.cent()

    def diffind(self,n,m,dra=10*pi/180,dta=10,nrot=41,nt=21,sr=0,stx=0,sty=0):
        """
        find most likely position of rod between images n and m
        dra,nrot: range and number of angles
        dta,nt: range and number of translate
        sr: offset angle of rotation
        stx,sty: offset of translation
        """
        rots = sr+linspace(-dra,dra,nrot)
        tt = linspace(-dta,dta,nt)
        r = array(self.rr,dtype=int)
        c = array(self.rc,dtype=int)

        lss = []
        lsr = []
        lsx = []
        lsy = []
        for rot in rots:
            corot = cos(rot)
            sirot = sin(rot)
            for tx in stx+tt:
                for ty in sty+tt:
                    dcc = self.dcc*corot-self.dcr*sirot+tx
                    dcr = self.dcc*sirot+self.dcr*corot+ty        
                    rr = array(self.cr + dcr,dtype=int)
                    rc = array(self.cc + dcc,dtype=int)
                    #print("rr = {}".format(rr))
                    #print("rc = {}".format(rc))                    
                    lss.append(mean((self.ims[n][r,c]-self.ims[m][rr,rc])**2))
                    lsx.append(tx)
                    lsy.append(ty)
                    lsr.append(rot)
        mini = lss.index(min(lss))
        return [lsr[mini]-sr,lsx[mini]-stx,lsy[mini]-sty]

    def seqapath(self,start=0,stop=-1,dra=10*pi/180,dta=10,nrot=41,nt=21):
        """tracking of original selected region"""
        self.rs = [0.0]
        self.xp = [self.cc]
        self.yp = [self.cr]
        
        if stop < 0:
            lind = len(self.ims)+stop
        else:
            lind = stop-1
        
        print("finding path in images (absolute)...")
        for k in range(start,lind):
            print ("step " + str(k))
            r,dx,dy = self.diffind(0,k+1,dra=dra,dta=dta,nrot=nrot,nt=nt,sr=self.rs[-1],stx=self.xp[-1]-self.xp[0],sty=self.yp[-1]-self.yp[0])
            self.rs.append(self.rs[-1]+r)
            self.xp.append(self.xp[-1]+dx)
            self.yp.append(self.yp[-1]+dy)
            print(str(k) + " to " + str(k+1))

    def seqpath(self,start=0,stop=-1,dra=10*pi/180,dta=10,nrot=41,nt=21):
        """tracking with update in each step - prone to dead-reckoning errors but able to adjust to change in appearance"""
        self.rs = [0.0]
        self.xp = [self.cc]
        self.yp = [self.cr]
        
        if stop < 0:
            lind = len(self.ims)+stop
        else:
            lind = stop-1
        
        print("finding path in images (relative)...")
        for k in range(start,lind):
            r,dx,dy = self.diffind(k,k+1,dra=dra,dta=dta,nrot=nrot,nt=nt)
            self.rs.append(self.rs[-1]+r)
            self.xp.append(self.xp[-1]+dx)
            self.yp.append(self.yp[-1]+dy)
            self.rodtransform(r,[dx,dy])
            print(str(k) + " to " + str(k+1))
        
    def showpath(self,n=-1,rectonly=False):
        im = copy(self.ims[n])
        if n < 0:
            n = len(self.rs)+n
        for k in range(n+1):
            if k>1 and not rectonly:
                lr,lc = line(int(self.yp[k-1]),int(self.xp[k-1]),int(self.yp[k]),int(self.xp[k]))
                im[lr,lc] = 1
            if k == n:
                sec = self.sec - self.xp[0]
                ser = self.ser - self.yp[0]
                corot = cos(self.rs[k])
                sirot = sin(self.rs[k])
                sec,ser = sec*corot-ser*sirot+self.xp[k],sec*sirot+ser*corot+self.yp[k]
                spr,spc = polygon_perimeter(ser,sec)
                im[spr,spc] = 1
        return im
        