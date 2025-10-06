import numpy as np
import functools
from math import sqrt, pi

def cmp(a, b):
    return float(a > b) - float(a < b)
    
class FccLattice:
    "Handles fcc crystal. Could be generalized to arbitrary crystal."
    def __init__(self, LatConst):
        self.a0 = np.array([0.5*LatConst,0.5*LatConst,0])
        self.a1 = np.array([0.5*LatConst,0,0.5*LatConst])
        self.a2 = np.array([0,0.5*LatConst,0.5*LatConst])
        Vc = np.dot(np.cross(self.a0,self.a1),self.a2) # Volume
        self.Volume = abs(Vc)
        print("Volume is", self.Volume) #print "Volume is ", self.Volume
        self.b0 = (2*pi/Vc)*np.cross(self.a1,self.a2)
        self.b1 = (2*pi/Vc)*np.cross(self.a2,self.a0)
        self.b2 = (2*pi/Vc)*np.cross(self.a0,self.a1)
        # Special points in Brillouin zone
        brs = 2*pi/LatConst
        self.GPoint = [0,0,0]
        self.LPoint = np.array([0.5,  0.5,  0.5])*brs
        self.KPoint = np.array([0.75, 0.75, 0])*brs
        self.XPoint = np.array([1.0,  0.0,  0])*brs
        self.WPoint = np.array([1.0,  0.5,  0])*brs
        
    def RMuffinTin(self):
        return 0.5*sqrt(np.dot(self.a0,self.a0)) # Spheres just touch

    def GenerateReciprocalVectors(self, q, CutOffK):
        # Many reciprocal vectors are generated and later only the shortest are used
        Kmesh0=[]
        for n in range(-q,q+1):
            for l in range(-q,q+1):
                for m in range(-q, q+1):
                    vec = n*self.b0+l*self.b1+m*self.b2
                    if np.dot(vec,vec) <= CutOffK**2:
                        Kmesh0.append(vec)
                    
        Kmesh0.sort(key=functools.cmp_to_key(lambda x, y: cmp(np.dot(x, x), np.dot(y, y))))
        self.Km = np.array(Kmesh0)
        print("K-mesh size =", len(self.Km))
        
    def ChoosePointsInFBZ(self, nkp, type=0): # Chooses the path in the 1BZ we will use
        
        def kv0(iq, q):
            return (iq-int((q+1.5)/2)+1)/(q+0.0)
        
        if type==0: # Choose mesh in the 1BZ to cover the whole space - for SC calculation
            kp=[]
            for i0 in range(nkp):
                r0 = kv0(i0,nkp)
                print('r0 =', r0)
                for i1 in range(nkp):
                    r1 = kv0(i1,nkp)
                    for i2 in range(nkp):
                        r2 = kv0(i2,nkp)
                        k = self.b0*r0+self.b1*r1+self.b2*r2
                        kp.append(k)
            print("Number of all k-points =", len(kp))

            kpc = []
            for k in kp:
                kpc.append(np.sort(k))
                
            # ChooseIrreducible k-points only
            irkp = []
            wkp  = []
            while len(kpc)>0: 
                tk = kpc[0]
                irkp.append(tk)
                wkp.append(0)
                # 48 symmetry operations of cubic group
                for ix in [-1,1]:
                    for iy in [-1,1]:
                        for iz in [-1,1]:
                            nk = np.sort([ix*tk[0], iy*tk[1], iz*tk[2]])
                            ii=0
                            while ii<len(kpc):
                                diff = np.sum(np.abs(nk - kpc[ii]))
                                if diff<1e-6:
                                    del kpc[ii]
                                    wkp[-1] += 1.
                                else:
                                    ii+=1

            self.wkp = np.array(wkp)/np.sum(wkp)
            self.kp = np.array(irkp)

            print("Number of irreducible k points is", len(self.kp))
            
        else:        # Choose one particular path in the 1BZ - for plotting purposes
            nkp = 4*int(nkp/4.)+1
            print("nkp =", nkp)
            self.kp = np.zeros((nkp,3), dtype=float)
            N0 = nkp // 4

            self.Points = [('$\\Gamma$', 0), ('$X$', N0), ('$L$', 2*N0), ('$\\Gamma$', 3*N0), ('$K$', 4*N0)]
            for i in range(N0): self.kp[i,:]      = self.GPoint + (self.XPoint-self.GPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0+i,:]   = self.XPoint + (self.LPoint-self.XPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*2+i,:] = self.LPoint + (self.GPoint-self.LPoint)*i/(N0-0.)
            for i in range(N0): self.kp[N0*3+i,:] = self.GPoint + (self.KPoint-self.GPoint)*i/(N0-0.)
            self.kp[4*N0] = self.KPoint
