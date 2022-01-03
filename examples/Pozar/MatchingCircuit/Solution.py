# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 11:01:35 2022

@author: hiro7
"""
"""
Purpose:
    Impedance matching and Tuning   
     Maximum power is delivered when the load is matched to the line (assuming the generator
is matched), and power loss in the feed line is minimized.
 Impedance matching sensitive receiver components (antenna, low-noise amplifier, etc.)
may improve the signal-to-noise ratio of the system.
 Impedance matching in a power distribution network (such as an antenna array feed
network) may reduce amplitude and phase errors.
"""
#input source impedance and load impedance, degree of freedom 
#Matching with lumpled element using L-network using two reactive element
#(Impedance-based/admittance based)
#whether to be used depends on the location of smith chart, namely if it's inside a circle of 1+jx,
#then impedance-based matching network is used

import numpy as np
import scipy as sp
import skrf as rf
import matplotlib.pyplot as plt
from matplotlib import cm
pi,exp,sin,cos,tan,sinh,cosh,tanh,sinc = [np.pi,np.exp,np.sin,np.cos,np.tan,\
                                     np.sinh,np.cosh,np.tanh,np.sinc]
inv,det = [np.linalg.inv,np.linalg.det]

sqrt = np.sqrt
arctan = np.arctan


def Matching_Zm(Zl,Z0): #the load impedance is inside the 1+1j, i.e. Rl is larger than Z0
    Rl,Xl = [Zl.real,Zl.imag]
    B1 = (Xl + sqrt(Rl/Z0)*sqrt(Rl**2 + Xl**2 - Z0*Rl))/(Rl**2 + Xl**2)
    B2 = (Xl - sqrt(Rl/Z0)*sqrt(Rl**2 + Xl**2 - Z0*Rl))/(Rl**2 + Xl**2)

    X1 = 1/B1 + Xl*Z0/Rl - Z0/B1/Rl
    X2 = 1/B2 + Xl*Z0/Rl - Z0/B2/Rl
    Zm1 = [X1,B1]
    Zm2 = [X2,B2]

    return [Zm1,Zm2]


def Matching_Ym(Zl,Z0): #the load impedance is outside the 1+1j, i.e. Rl is smaller than Z0
    Rl,Xl = [Zl.real,Zl.imag]
    X1 = sqrt(Rl*(Z0-Rl)) -Xl
    X2 = -sqrt(Rl*(Z0-Rl)) - Xl
    
    B1 = sqrt((Z0 - Rl)/Rl)/Z0
    B2 = -sqrt((Z0 - Rl)/Rl)/Z0
    Ym1 = [X1,B1]
    Ym2 = [X2,B2]

    return [Ym1,Ym2]


    
"""
Shunt stubs are
preferred for microstrip line or stripline, while series stubs are preferred for slotline or
coplanar waveguide.
two adjustable parameters are the distance, d, from the load
to the stub position, and the value of susceptance or reactance provided by the stub.
"""


test1 = Matching(200-100j,100)

"""
Both combinations are available, but:
    One solution, however, may
result in significantly smaller values for the reactive components, or may be the preferred
solution if the bandwidth of the match is better, or if the SWR on the line between the
matching network and the load is smaller.
"""

from skrf.media import MLine
freq = rf.Frequency(2,2,1,unit = 'ghz') #2 GHz
Line = MLine(frequency = freq, w = 0.483e-3,h = 0.5e-3,t = 0.01e-3, ep_r = 9.9) #Alumina line

def shunt_stub(Z0,Zl,Frequency,Line,stub_type = 'open'): #Y0 is characteristic impedance of stub
    c0 = 2.997924e8
    Rl,Xl = [Zl.real,Zl.imag]    
    #Shunt Stub
    t1 = (Xl + sqrt(Rl/Z0*((Z0 - Rl)**2 + Xl**2)))/(Rl - Z0) 
    t2 = (Xl - sqrt(Rl/Z0*((Z0 - Rl)**2 + Xl**2)))/(Rl - Z0) 
    t = np.array([t1,t2])
    d = arctan(t)/Line.beta #d is a distance and beta can be derived by line parameter
    
    #G = Rl*(1 + t**2)/(Rl**2 + (Xl + Z0*t)**2)
    B = (Rl**2*t - (Z0 - Xl*t)*(Xl + Z0*t))/Z0/(Rl**2 + (Xl + Z0*t)**2)
    lambda0 = c0/freq.f
    if stub_type == 'open':
        l = -lambda0/(2*pi)*arctan(Line.z0*B)
    else: #short
        l = -lambda0/(2*pi)*arctan(1/Line.z0/B)
        
    return [d,l]
