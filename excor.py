""" This module defines class to compute
exchange-correlation potential and
exchange-correlation energy within LDA

class name:
    ExchangeCorrelation(type)  --     type can be 1,2,3,4,5
Adopted from Kristjan Haule
Modified more and optimized by Edmond Febrinicko Armay    
"""
from math import pi
from math import exp
from math import sin
from math import sqrt
from math import atan
from math import log
from pythtb import *
from numpy import array
#from math import cos
class  ExchangeCorrelation(object):
    """******************************************************************************/
    Calculates Exchange-Correlation Energy and Potential                       */ 
    type=0 - due to U.von.Barth and L.Hedin, J.Phys.C5, 1629 (1972)            */
    type=1 - O.E.Gunnarsson and S.Lundqvist,  Phys.Rev.B                       */
    type=2 - V.L.Moruzzi, J.F.Janak, and A.R.Williams, Calculated              */
             Electronic Properties of Metals (New York, Pergamon Press, 1978)  */
    type=3 - S.H.Vosko, L.Wilk, and M.Nusair, Can.J.Phys.58, 1200 (1980)       */
    type=4 - Correlation of Perdew and Wang 1991                               */
    type=5 - Correlation of Chachiyo                                           */
    ******************************************************************************/
    """
    #def __init__(self, type_=3):
    def __init__(self, type_=5):
        self.type = type_
        #self.alphax = 0.610887057710857 #//(3/(2 Pi))^(2/3)
        self.alphax = pow(3/(2*pi), 2/3)
        self.Aw = 0.0311 
        self.Bw = -0.048 
        self.Cw = 0.002 
        self.D  = -0.0116 
        self.gamma  = -0.1423 
        self.beta1  =  1.0529 
        self.beta2  =  0.3334 
        self.Ap  =  0.0621814 
        self.xp0 = -0.10498 
        self.bp  =  3.72744 
        self.cp  =  12.9352 
        self.Qp  =  6.1519908 
        self.cp1 =  1.2117833 
        self.cp2 =  1.1435257 
        self.cp3 = -0.031167608
        #self.kF = 1
        if self.type==0:
            self.C = 0.0504;
            self.A = 30
        if self.type==1:
            self.C = 0.0666;
            self.A = 11.4
        if self.type==2:
            self.C = 0.045; 
            self.A = 21

    def Vx(self, rs): # Vx
        copper=w90(r"copper")
        fermi_ev=0.7E+01
        my_model=copper.model(zero_energy=fermi_ev, min_hopping_norm=0.01, max_distance=None, ignorable_imaginary_part=0.01)
        
        #return -self.alphax/rs
        #return (-1/(2*pow(pi, 2)))*(exp(-self.alphax*rs)/pow(rs, 3))*((sin(self.alphax*rs)/rs)-(self.alphax*cos(self.alphax*rs)))
        return (-1/(2*pow(pi, 2)))*(exp(-self.alphax*rs)/pow(rs, 3))*(sin(self.alphax*rs)/rs) #Yukawa's screened exchange
        
    def ExVx(self, rs): # Ex-Vx
        return 0.25*self.alphax/rs
    
    """def Ex(self, rs):
        return -0.75*self.alphax/rs"""

    def Vc(self, rs): # Vc
        if (self.type<3):
            x = rs/self.A;
            return -0.5*self.C*log(1+1./x)
        elif(self.type<4): # type=3 VWN
            x=sqrt(rs)
            xpx= x*x + self.bp*x + self.cp
            atnp = atan(self.Qp/(2*x+self.bp))
            ecp = 0.5*self.Ap*(log(x*x/xpx)+self.cp1*atnp-self.cp3*(log((x-self.xp0)**2/xpx)+self.cp2*atnp))
            return ecp - self.Ap/6.*(self.cp*(x-self.xp0)-self.bp*x*self.xp0)/((x-self.xp0)*xpx)
        #else:
        elif(self.type<5):
            if rs>1:
                return self.gamma/(1+self.beta1*sqrt(rs)+self.beta2*rs)*(1+7/6.*self.beta1*sqrt(rs)+self.beta2*rs)/(1+self.beta1*sqrt(rs)+self.beta2*rs)
            else:
                return self.Aw*log(rs)+self.Bw-self.Aw/3.+2/3.*self.Cw*rs*log(rs)+(2*self.D-self.Cw)*rs/3.
        else:
            #paramagnetic
            a=(log(2)-1)/(2*pow(pi, 2))
            b=20.4562557
            #ferromagnetic
            #a=(log(2)-1)/(4*pow(pi, 2))
            #b=27.4203609
            return a*log(1+(b/rs)+(b/pow(rs, 2)))
    def EcVc(self, rs): # Ec-Vc
        if self.type<3 :
            x = rs/self.A
            epsilon = -0.5*self.C*((1+x*x*x)*log(1+1/x)+0.5*x-x*x-1/3.)
            return epsilon-Vc(rs)
        elif self.type<4: # type=3 VWN
            x = sqrt(rs)
            return self.Ap/6.*(self.cp*(x-self.xp0)-self.bp*x*self.xp0)/((x-self.xp0)*(x*x+self.bp*x+self.cp))
        #else:
        elif self.type<5:
            if rs>1:
                return 2*self.gamma/(1+self.beta1*sqrt(rs)+self.beta2*rs)-Vc(rs)
            else:
                return self.Aw*log(rs)+self.Bw+self.Cw*rs*log(rs)+self.D*rs-Vc(rs)
        else:
            return 0