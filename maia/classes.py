import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from skrf.media import DefinedGammaZ0, Freespace
import schemdraw
import schemdraw.elements as elm


from .functions import *

pi = np.pi
sqrt = np.sqrt
sin = np.sin
cos = np.cos

jv = sp.special.jv
kv = sp.special.kv
kn = sp.special.kn
iv = sp.special.iv
jvp = sp.special.jvp
kvp = sp.special.kvp
hankel2 = sp.special.hankel2
fminbound = sp.optimize.fminbound #This is used to determine u in dispersion relations to find x-value giving minimum y
gradient = np.gradient
e0 = 8.85e-12
m0 = 4*pi*1e-7
c0 = 2.9979e8
atan = np.arctan



"""
Filter Analysis
from the specification such as center frequency, fractional bandwidth, etc...
Prototype filter with g-parameters are assigned
def Filter(g,fc,filter_type = 'L_base'):
def BPF(g,fc,FBW,filter_type = 'L_base'):
    
output Network class of Filter and LC value

"""


#g_butterworth(N) and output is a list

class LowPassFilter(object):
    """
    test
    """
    def __init__(self,N = 2, fc = 1e9, Base = 'L_base', Type = 'Butterworth',ripple = None,f_start = 1,f_stop=10,Npoints=1001,Zl = 50):
        self.f_start = f_start
        self.f_stop = f_stop
        self.Npoints = Npoints
        self.Zl = Zl
        self.ripple = ripple
        self.Freq = rf.Frequency(self.f_start,self.f_stop,self.Npoints,'ghz')
        self.N = N
        self.fc = fc
        self.Type = Type 
        self.Base = Base
        if Type == 'Butterworth':
            self.g = g_butterworth(self.N)[1:-1]
            self.load = g_butterworth(self.N)[-1]*self.Zl
        elif Type == 'Chebyshev':
            self.g = g_chebyshev(self.N,self.ripple)[1:-1]
            self.load = g_chebyshev(self.N,self.ripple)[-1]*self.Zl
        self.Network,self.LC = LPF(self.g,self.fc,filter_type=self.Base,Freq = self.Freq)
        self.L = self.LC[0]
        self.C = self.LC[1]
        
    #added methods: plot_summary(),Implemantation() using  inverters...            
    def plot_summary(self,f_start = None,f_stop = None, ax = ((None,None),(None,None))):
        if ax == ((None,None),(None,None)):
            fig,ax = plt.subplots(2,2,figsize = (8,8))
        if (f_start,f_stop) == (None,None):
            f_start = self.f_start
            f_stop = self.f_stop
        freq_range = str(f_start)+'-'+str(f_stop)+' GHz'
        
        self.Network.s21[freq_range].plot_s_db(ax=ax[0,0],lw=2)
        self.Network.s11[freq_range].plot_s_db(ax=ax[0,1],lw=2)
        
        ax[0,0].set_ylabel('Loss [dB]')
        ax[0,1].set_ylabel('Return Loss [dB]')
        #ax1.set_ylim(-1,0)
        #ax1_2.set_ylim(-30,0)
        #ax[0][0].legend().remove()
        #ax[0][1].legend().remove()
        
        self.Network.s11[freq_range].plot_s_smith(ax = ax[1,0],lw = 2)
        ax[1,1].plot(self.Network.frequency[freq_range].f,\
                     self.Network.s21[freq_range].group_delay[:,0,0]*1e9)
        #ax[1].set_xlabel('Frequency [GHz]')        
        #if ax[1].get_ylabel != 'Group Delay [ns]': 
        #    ax[1].set_ylabel('Group Delay [ns]')

        
        
    def StepImpedance(self,Z,h,t,ep_r,tanD):
        S,length,Line,w,beta = StepImpedanceLPF(self.N,self.g,Z,h,t,ep_r,tanD,self.Freq,self.fc/1e9)
        d = Step_Schematic(Z,Step_test[1])
        return S,length,Line,w,beta,d 
   
    
   
    def schematic(self):
      d = Schematic_LPF(self.LC,self.Base)
      return d  




        
        

        
    """    
    def implementation(self,Z0,[h,t,ep_r,tanD]):
        StepImpedanceLPF(self.N,self.g,Z0,h,t,ep_r,tanD) #Z0 = [Zsource,Zhigh,Zlow]
    """    
        
class BandPassFilter(LowPassFilter):
    def __init__(self,N = 5, fc = 1e9, FBW = 0.1, Base = 'L_base', Type = 'Butterworth',ripple = None,f_start = 1,f_stop=10,Npoints=1001,Zl = 50):
        super().__init__(N,fc,Base,Type,ripple,f_start,f_stop,Npoints,Zl)
        self.FBW = FBW
        self.Network,self.LC = BPF(self.g,self.fc,self.FBW,filter_type=self.Base,Freq = self.Freq)
        self.L = self.LC[0]
        self.C = self.LC[1]
 
    def schematic(self):
        d = Schematic_BPF(self.LC,self.Base)
        return d
    
    #def plot_summary(self,f_start = None,f_stop = None, ax = ((None,None),(None,None))):
    #    super()

    """
    def plot_summary(self,f_start = None,f_stop = None):
        if [f_start,f_stop] == [None,None]:
            f_start = self.f_start
            f_stop = self.f_stop
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1_2 = ax1.twinx()
        
        self.Network.s21.plot_s_db(ax=ax1,color='orange',lw=2)
        self.Network.s11.plot_s_db(ax=ax1_2,color='blue',lw=2)
        
        ax1.set_ylabel('Loss [dB]')
        ax1_2.set_ylabel('Return Loss [dB]')
        ax1.set_ylim(-1,0)
        ax1_2.set_ylim(-30,0)
        ax1.set_xlim(2e9,3e9)
        ax1_2.set_xlim(2e9,3e9)
        ax1.legend().remove()
        ax1_2.legend().remove()
        ax2 = fig.add_subplot(212)
        ax2.plot(self.Network.f/1e9,self.Network.s21.group_delay[:,0,0]*1e9,color='orange',lw = 3)
        ax2.set_xlim(2,3)
    
        ax2.set_xlabel('Frequency [GHz]')        
        ax2.set_ylabel('Group Delay [ns]')
    """        
"""
class Test0():
    def __init__(self,a,b):
        self.a = a
        self.b = b
    def add(self):
        return self.a + self.b
        
class Test1(Test0):
    def __init__(self,a,b,c):
        super().__init__(a,b)
        self.c = c
"""

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
        .solve_reactance()
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

    def __init__(self,Z0,Zl,fc,fbw,npoints):
        c0 = 2.99792458e8
        f2beta0 = lambda f:2*pi*f/c0
        # load Network
        self.Z0 = Z0
        self.Frequency = rf.Frequency((fc-fc*fbw/2)/1e9,(fc+fc*fbw/2)/1e9,npoints,unit = 'ghz')
        self.fc = self.Frequency.center
        self.media = DefinedGammaZ0(frequency = self.Frequency, z0= self.Z0,gamma = 1j*f2beta0(self.Frequency.f))
        self.Zl = Zl
        self.load = self.media.load(rf.zl_2_Gamma0(self.Z0,self.Zl))
        #media.load(rf.zl_2_Gamma0(100,200-100j)).plot_s_db()
        self.Rl,self.Xl = [Zl.real,Zl.imag] #lumped element
        self.XB1,self.XB2,self.circuit_type = self.solve_reactance() #reactance and susceptance
        self.LC1,self.LC2 = [self.find_lumped(self.XB1),self.find_lumped(self.XB2)]
        #LC value, flag, Networkcalculated inductance and capacitance to output of Network file
        self.schematic1 = self.schematic_LC(self.LC1[0],self.LC1[1])
        self.schematic2 = self.schematic_LC(self.LC2[0],self.LC2[1])
        self.schematic_size = self.schematic1.get_bbox()
        self.Network1,self.Network2 = [self.LC1[-1],self.LC2[-1]]
        self.matched1,self.matched2 = [self.Network1**self.load,self.Network2**self.load]
    def solve_reactance(self): #LC or CL
        
        if self.Rl >= self.Z0:
            B1 = (self.Xl + sqrt(self.Rl/self.Z0)*sqrt(self.Rl**2 + self.Xl**2 - self.Z0*self.Rl))/(self.Rl**2 + self.Xl**2)
            B2 = (self.Xl - sqrt(self.Rl/self.Z0)*sqrt(self.Rl**2 + self.Xl**2 - self.Z0*self.Rl))/(self.Rl**2 + self.Xl**2)
        
            X1 = 1/B1 + self.Xl*self.Z0/self.Rl - self.Z0/B1/self.Rl
            X2 = 1/B2 + self.Xl*self.Z0/self.Rl - self.Z0/B2/self.Rl
            #self.Zm1 = [self.X1,self.B1]
            #self.Zm2 = [self.X2,self.B2]
            circuit_type = 'impedance'

            
            
        else: #circuit_type == 'Admittance'
            X1 = sqrt(self.Rl*(self.Z0-self.Rl)) - self.Xl
            X2 = -sqrt(self.Rl*(self.Z0-self.Rl)) - self.Xl
            
            B1 = sqrt((self.Z0 - self.Rl)/self.Rl)/self.Z0
            B2 = -sqrt((self.Z0 - self.Rl)/self.Rl)/self.Z0
            circuit_type = 'admittance'        
            
        return [X1,B1],[X2,B2],circuit_type
        
        
        
    def find_lumped(self,XB): #[X, B]
        if self.circuit_type == 'impedance': #from Zm to LC
            if XB[0] > 0: #X is inductance otherwise, capacitance                
                LC = [XB[0]/(2*pi*self.fc),XB[1]/(2*pi*self.fc)]#LC
                LCflag = True
                Network = self.media.inductor(LC[0]) ** self.media.shunt_capacitor(LC[1])
                
            else:
                LC = [-1/(XB[0]*2*pi*self.fc),-1/(XB[1]*2*pi*self.fc)]#CL
                LCflag = False
                Network = self.media.capacitor(LC[0])**self.media.shunt_inductor(LC[1])
            
        if self.circuit_type == 'admittance': #from Ym to LC
            if XB[0] > 0: #X is inductance otherwise, capacitance                
                LC = [XB[0]/(2*pi*self.fc),XB[1]/(2*pi*self.fc)]#CL
                LCflag = True
                Network = self.media.shunt_capacitor(LC[0])**self.media.inductor(LC[1])
            else:
                LC = [-1/(XB[0]*2*pi*self.fc),-1/(XB[1]*2*pi*self.fc)]#LC
                LCflag = False
                Network =self.media.shunt_inductor(LC[0])**self.media.capacitor(LC[1])
        return LC,LCflag,Network
                
                
    def schematic_LC(self,LC,flag):
        if self.circuit_type == 'impedance':
            d = schemdraw.Drawing()    
            Start = d.add(elm.Resistor().label(str(self.Z0) + ' Ohm'))   
            if flag == True:
                d += elm.Inductor().label(str(np.round(self.unit_change(LC[0])[0],2)) + self.unit_change(LC[0])[1] +'H') 
                d.push()
                d += elm.Capacitor().down().label(str(np.round(self.unit_change(LC[1])[0],2)) + self.unit_change(LC[1])[1] +'F',rotate = True,loc = 'bottom') 
                d.pop()
                d += elm.Line()
                d += elm.Resistor().down().label(str(self.Zl) + ' Ohm',rotate = True,loc = 'bottom')
                d += elm.Line().left().tox(Start.start)
            else:
                d += elm.Capacitor().label(str(np.round(self.unit_change(LC[0])[0],2)) + self.unit_change(LC[0])[1] +'F') 
                d.push()
                d += elm.Inductor().down().label(str(np.round(self.unit_change(LC[1])[0],2)) + self.unit_change(LC[1])[1] +'H',rotate = True,loc = 'bottom') 
                d.pop()
                d += elm.Line()
                d += elm.Resistor().down().label(str(self.Zl) + ' Ohm',rotate = True,loc = 'bottom')
                d += elm.Line().left().tox(Start.start)            
            

        elif self.circuit_type == 'admittance':
            d = schemdraw.Drawing()    
            Start = d.add(elm.Resistor().label(str(self.Z0) + ' Ohm'))   
            if flag == True:
                d.push()
                d += elm.Capacitor().down().label(str(np.round(self.unit_change(LC[0])[0],2)) + self.unit_change(LC[0])[1] +'F',rotate = True,loc = 'bottom') 
                d.pop()
                d += elm.Inductor().label(str(np.round(self.unit_change(LC[1])[0],2)) + self.unit_change(LC[1])[1] +'H') 
                d += elm.Resistor().down().label(str(self.Zl) + ' Ohm',rotate = True,loc = 'bottom')
                d += elm.Line().left().tox(Start.start)
            else:
                d.push()
                d += elm.Inductor().down().label(str(np.round(self.unit_change(LC[0])[0],2)) + self.unit_change(LC[0])[1] +'H',rotate = True,loc = 'bottom') 
                d.pop()
                d += elm.Capacitor().label(str(np.round(self.unit_change(LC[1])[0],2)) + self.unit_change(LC[1])[1] +'F') 
                d += elm.Resistor().down().label(str(self.Zl) + ' Ohm',rotate = True,loc = 'bottom')
                d += elm.Line().left().tox(Start.start)
        return d
    
    def plot_summary(self,figsize = (8,4),hspace = 0.2,wspace = 0):
        
        fig,ax = plt.subplots(1,2,figsize = figsize)
        fig2,ax2 = plt.subplots(2,1,figsize = (6.4,5.4))
        fig.subplots_adjust(hspace = hspace,wspace = wspace)
        ax2[0].axis('off')
        ax2[1].axis('off')
        xmin,ymin,xmax,ymax = self.schematic_size
        self.matched1.plot_s_mag(ax = ax[0],label = 'Sol1') 
        self.matched2.plot_s_mag(ax = ax[0],label = 'Sol2')
        self.matched1.plot_s_smith(ax = ax[1],label = 'Sol1') 
        self.matched2.plot_s_smith(ax = ax[1],label = 'Sol2')

        self.schematic1.draw(ax = ax2[0])
        self.schematic2.draw(ax = ax2[1])
        ax2[0].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax2[1].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        fig2.subplots_adjust(hspace = 0.2)
        

    
    def unit_change(self,A):
        unit = None
        if A >= 1e-18 and A < 1e-15:
            A = A*1e18
            unit = 'a'
        elif A >= 1e-15 and A < 1e-12:
            A = A*1e15
            unit = 'f'
        elif A >= 1e-12 and A < 1e-9:
            A = A*1e12
            unit = 'p'
        elif A >= 1e-9 and A < 1e-6:
            A = A*1e9
            unit = 'n'
        elif A >= 1e-6 and A < 1e-3:
            A = A*1e6
            unit = 'u'
        return A,unit
        
    def define_media(self,w,h,t,Er,tanD): #microstril line or other line defined here
        self.Line = MLine(frequency = self.Frequency, w = w,h = h,t = t, ep_r = Er,tand = tanD,z0 = self.Z0) #Alumina line

    def shunt_stub(self,stub_type = 'open'): #Y0 is characteristic impedance of stub
        #Shunt Stub
        t1 = (self.Xl + sqrt(self.Rl/self.Z0\
                             *((self.Z0 - self.Rl)**2 + self.Xl**2)))/(self.Rl - self.Z0) 
        t2 = (self.Xl - sqrt(self.Rl/self.Z0\
                             *((self.Z0 - self.Rl)**2 + self.Xl**2)))/(self.Rl - self.Z0) 
        #t = np.array([t1,t2])
        d1 = arctan(t1)/(2*pi*self.fc)*c0 #d is a distance and beta can be derived by line parameter
        d2 = arctan(t2)/(2*pi*self.fc)*c0 #d is a distance and beta can be derived by line parameter
        line_extend1 = self.Line.line(d1,unit = 'm')
        line_extend2 = self.Line.line(d2,unit = 'm')
        #G = Rl*(1 + t**2)/(Rl**2 + (Xl + Z0*t)**2)
        B1 = (self.Rl**2*t1 - (self.Z0 - self.Xl*t1)*(self.Xl + self.Z0*t1))/self.Z0/(self.Rl**2 + (self.Xl + self.Z0*t1)**2)
        B2 = (self.Rl**2*t2 - (self.Z0 - self.Xl*t2)*(self.Xl + self.Z0*t1))/self.Z0/(self.Rl**2 + (self.Xl + self.Z0*t2)**2)

        lambda0 = c0/self.fc
        if stub_type == 'open':
            l1 = -lambda0/(2*pi)*arctan(self.Line.z0*B1)
            l2 = -lambda0/(2*pi)*arctan(self.Line.z0*B2)

            stub_shunt = [self.Line.shunt_delay_open(i,unit = 'm') for i in [l1,l2]]
        else: #short
            l1 = -lambda0/(2*pi)*arctan(1/self.Line.z0/B1)
            l2 = -lambda0/(2*pi)*arctan(1/self.Line.z0/B2)
            stub_shunt = [self.Line.shunt_delay_short(i,unit = 'm') for i in [l1,l2]]

        #shunt_line
        stub1 = stub_shunt[0]**line_extend1
        stub2 = stub_shunt[1]**line_extend2
        
        
        self.stub_matched1 = stub1 ** self.load
        self.stub_matched2 = stub2 ** self.load
        #self.matched1,self.matched2 = [self.Network1**self.load,self.Network2**self.load]

        return [self.stub_matched1,self.stub_matched2]
    #def quarter(self): #single section of quarter wavelength
    """
    One drawback of the quarter-wave transformer is that it can only match a real load
    impedance. A complex load impedance can always be transformed into a real impedance,
    however, by using an appropriate length of transmission line between the load and the
    transformer, or an appropriate series or shunt reactive element.
    """
    #    return 0
    

        #def schematic_lumped(self):
    """
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
    """  
            
