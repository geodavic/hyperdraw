import numpy as np
import math
import cmath
import sys

#This script meant to do the same as mathematica script of same name, but never tested this one
#requires: python3, numpy, 


#octagon class
class polygon:
    def __init__(self,pts,labels=None):
        self.n = len(pts)
        self.pts = [HNormalize(p) for p in pts]
        self.lins = np.array([np.array([pts[i],pts[(i+1)%self.n]]) for i in range(len(pts))])
        if labels is None: 
            self.labels = ["a","b","a","b","c","d","c","d"]
        else:
            self.labels = labels

    def reflect(self):
        octs =[]
        for l in self.lins:     
            newpts = np.array([reflect(p,line(l[0],l[1])) for p in self.pts])
            O = polygon(newpts,labels=self.labels)
            octs.append(O)
        return octs

#stereographic projection to Poincar\'{e} disk
def pip(p):
    c = 1.0*(p[2]+1)
    return np.array([p[0],p[1]])/c

#inverse stereographic projection
def pipInv(v):
    a,b = v[0],v[1]
    c = 1.0*(1-a**2-b**2)
    return np.array([2*a,2*b,(1+a**2+b**2)])/c

#stereographic projection from dual space to complement of Poincar\'{e} disk
def pil(l):
    c = 1.0*(-l[2])
    return np.array([l[0],l[1]])/c

#normalize to hyperbola
def HNormalize(p):
    return p/np.sqrt(abs(-p[0]**2-p[1]**2+p[2]**2))

#line determined by pair of points (as element of dual space)
def line(v1,v2):
    l=np.cross(v1,v2)
    return HNormalize(l)
  
#reflect a point about a line
def reflect(p,L):
    dist = np.linalg.norm(pip(p)-pil(L))/abs(-1/L[2])
    p2 = pil(L) - (pil(L)-pip(p))/(dist**2)
    return pipInv(p2)

#False iff p is too close to the boundary of the disk
def filterPoint(p,thresh=0.99):
    return np.linalg.norm(pip(p)) < thresh

def draw(scale,recLim=2):
    sincos = math.sqrt(2)/2
    #rotation matrix by pi/4
    r = np.array([[sincos,-sincos,0],[sincos,sincos,0],[0,0,1]])
    pS = np.array([0,scale])
    p = pipInv(pS)
    pts = [p]
    for i in range(7):
        pts.append(r@p)
        p = pts[-1]
    pts = np.array(pts)
    #generate reflected octagons
    octagons = [polygon(pts)]
    for i in range(recLim):
        size = len(octagons)
        for j in range(size):
            octagons += octagons[j].reflect()

    #lines to draw
    lines2draw = []
    labels = []
    for i in range(len(octagons)):
        lines2draw += list(octagons[i].lins)
        labels += octagons[i].labels
    
    #remove lines that are too small
    lines2drawFiltered = []
    for l in lines2draw:
        if filterPoint(l[0]) or filterPoint(l[1]):
            sorted_l = l[l[:,0].argsort()]
            lines2drawFiltered += [np.around(sorted_l,decimals=10)] 
    
    lines2drawFiltered = np.array(lines2drawFiltered)
    uniquelines,indices = np.unique(lines2drawFiltered,axis=0,return_index=True)
    return uniquelines,[labels[idx] for idx in indices]

#export lines and labels for mathematica input
def export(lines,labels):
    lines.flatten().astype('float32').tofile('lines.dat')
    fp = open("./labels.csv","w")
    s = ""
    for l in labels:
        s+= "{},".format(l)
    s = s[:-1]
    fp.write(s)
    fp.close()

lins,labels = draw(float(sys.argv[1]),recLim=int(sys.argv[2]))
export(lins,labels)

