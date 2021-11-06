import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from skrf.media import DefinedGammaZ0, Freespace
import schemdraw
import schemdraw.elements as elm


from .maia_functions import *

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
    def __init__(self,N = 2, fc = 1e9, BW = None, Base = 'L_base', Type = 'Butterworth',ripple = None,f_start = 1,f_stop=10,Npoints=1001,Zl = 50):
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
    def plot_summary(self,f_start = None,f_stop = None, ax = [[None,None],None]):
        if ax == ((None,None),None):
            fig,ax = plt.subplots(2,1)
            ax1_2 = ax[0].twinx()
            ax = ((ax[0],ax1_2),ax[1])
        if [f_start,f_stop] == [None,None]:
            f_start = self.f_start
            f_stop = self.f_stop
        
        
        self.Network.s21[str(f_start)+'-'+str(f_stop)+' GHz'].plot_s_db(ax=ax[0][0],color='orange',lw=2)
        self.Network.s11[str(f_start)+'-'+str(f_stop)+' GHz'].plot_s_db(ax=ax[0][1],color='blue',lw=2)
        
        ax[0][0].set_ylabel('Loss [dB]')
        ax[0][1].set_ylabel('Return Loss [dB]')
        #ax1.set_ylim(-1,0)
        #ax1_2.set_ylim(-30,0)
        ax[0][0].legend().remove()
        ax[0][1].legend().remove()
        
        ax[1].plot(self.Network.frequency[str(f_start)+'-'+str(f_stop)+' GHz'].f,self.Network[str(f_start)+'-'+str(f_stop)+' GHz'].s21.group_delay[:,0,0]*1e9,color='orange',lw = 3)
    
        ax[1].set_xlabel('Frequency [GHz]')        
        if ax[1].get_ylabel != 'Group Delay [ns]': 
            ax[1].set_ylabel('Group Delay [ns]')

        
        
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
        
class BandPassFilter():
    def __init__(self,N = 5, fc = 1e9, FBW = 0.1, Base = 'L_base', Type = 'Butterworth',ripple = None,f_start = 1,f_stop=10,Npoints=1001,Zl = 50):
        self.N = N
        self.Npoints = Npoints
        self.f_start = f_start
        self.ripple = ripple
        self.f_stop = f_stop
        self.Freq = rf.Frequency(self.f_start,self.f_stop,self.Npoints,'ghz')
        self.fc = fc
        self.Type = Type 
        self.Zl = Zl
        self.Base = Base
        if Type == 'Butterworth':
            self.g = g_butterworth(self.N)[1:-1]
            self.load = g_butterworth(self.N)[-1]*self.Zl
        elif Type == 'Chebyshev':
            self.g = g_chebyshev(self.N,self.ripple)[1:-1]
            self.load = g_chebyshev(self.N,self.ripple)[-1]*self.Zl

        self.FBW = FBW
        self.Network,self.LC = BPF(self.g,self.fc,self.FBW,filter_type=self.Base,Freq = self.Freq)
        self.L = self.LC[0]
        self.C = self.LC[1]
        
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
        
    def schematic(self):
        d = Schematic_BPF(self.LC,self.Base)
        return d
        


