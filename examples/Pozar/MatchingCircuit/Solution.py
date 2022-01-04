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
class MatchingCircuit(object):
    """
    class MatchingCircuit(Z0,Zl,frequency)
    input parameters are following:
        Z0: source impedance, real value
        Zl: load impedance, complex value
        frequency: matching frequency
        
    method
        *initial methods
        .analytical()
            calculates analytical reactance X and susceptance B of the matching circuit
            there are two types of circuit, impedance circuit or admittance circuit
            in accordance with the relation between the load impedance and the source impedance.
            if the load impedance is larger than source impedance, i.e. located at outside 1+1j circle in a smith chart,
            impedance circuit is chosen. If not, admittance circuit.
            
        .lumped()
            determines lumped element value from analytical solutions above
            generates also schematic class using schemdraw
            
        *others
        .stub()
            with given line information such as characteristic impedance and propagation constant of microstrip line or waveguide,
            single-stub matching circuit is constructed
        .quarter()
            with given line information such as characteristic impedance and propagation constant of microstrip line or waveguide,
            a single quarter-wavelength circuit is constructed
        
    """

    def __init__(self,Z0,Zl,frequency):
        self.Z0 = Z0
        self.Rl,self.Xl = [Zl.real,Zl.imag] #lumped element
        self.XB1,self,XB2,self.circuit_type = self.analytical() #reactance and susceptance
        self.lumped_expansion()
        #calculated inductance and capacitance to output of Network file
        

    def analytical(self): #LC or CL
        if self.Rl >= self.Z0:
            self.B1 = (self.Xl + sqrt(self.Rl/self.Z0)*sqrt(self.Rl**2 + self.Xl**2 - self.Z0*self.Rl))/(self.Rl**2 + self.Xl**2)
            self.B2 = (self.Xl - sqrt(self.Rl/self.Z0)*sqrt(self.Rl**2 + self.Xl**2 - self.Z0*self.Rl))/(self.Rl**2 + self.Xl**2)
        
            self.X1 = 1/self.B1 + self.Xl*self.Z0/self.Rl - self.Z0/self.B1/self.Rl
            self.X2 = 1/self.B2 + self.Xl*self.Z0/self.Rl - self.Z0/self.B2/self.Rl
            #self.Zm1 = [self.X1,self.B1]
            #self.Zm2 = [self.X2,self.B2]
            self.circuit_type = 'impedance'
            return [self.X1,self.B1],[self.X2,self.B2],self.circuit_type
            
            
        else: #circuit_type == 'Admittance'
            self.X1 = sqrt(self.Rl*(self.Z0-self.Rl)) - self.Xl
            self.X2 = -sqrt(self.Rl*(self.Z0-self.Rl)) - self.Xl
            
            self.B1 = sqrt((self.Z0 - self.Rl)/self.Rl)/self.Z0
            self.B2 = -sqrt((self.Z0 - self.Rl)/self.Rl)/self.Z0
            self.circuit_type = 'admittance'        
            return [self.B1,self.X1],[self.B2,self.X2],self.circuit_type
        
        
        
        def lumped(self,XB):
            if self.circuit_type == 'impedance': #from Zm to LC
                L1,L_flag = self.X2LC(XB[0])
                C2,C_flag = self.B2LC(XB[1])
                
                if L_flag == True:
                    self.L = skrf.Inducgtor()
                    self.C = skrf.shunt_Capacitor()
                    self.LC = self.L**self.C
                else:
                    self.C = skrf.Capacitor()
                    self.L = skrf.shunt_Inductor
                    self.LC = self.C**self.L                    
                
            if self.circuit_type == 'admittance': #from Ym to LC
                C1,C_flag = self.B2LC(XB[0])
                L2,L_flag = self.X2LC(XB[1])     
                
                if C_flag == True:
                    self.C = skrf.shunt_capacitor()
                    self.L = skrf.inductor()
                    self.LC = self.C ** self.L
                else:
                    self.L = skrf.shunt_inductor()
                    self.C = skrf.capacitor()
                    self.LC = self.L ** self.C
                    
                    
        def schematic_lumped(self):
            
        def X2LC(self,X):
            if x > 0: #
                L_flag = True
                return X/(2*pi*self.frequency),L_flag
            else: #C
                L_flag= False
                return 1/(2*pi*self.frequency*X),L_flag
            
        def B2LC(self,B):
            if B > 0:
                C_flag = True
                return B/(2*pi*self.frequency),C_flag
            else:
                C_flag = False
                return 1/(2*pi*self.frequency*B),C_flag
            
        
    def stub(Line): #series or shunt
        
test = MatchingCircuit(100,50+20j,500e6)





    
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
