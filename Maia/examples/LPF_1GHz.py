import os
import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,brentq,minimize, differential_evolution
from skrf.media import DefinedGammaZ0,MLine,Freespace
import schemdraw
import schemdraw.elements as elm
path = r"C:\Users\hiro7\Documents\Python\Maia\Package"
os.chdir(path)

#Reference
#https://schemdraw.readthedocs.io/en/latest/usage/placement.html


from Maia_functions import *
import Maia

"""
Example1
Design a maximally flat low-pass filter with a cutoff frequency of 2 GHz, impedance
of 50, with stages N of 5, so that at least 15 dB insertion loss at 3 GHz.  compare with an
equal-ripple (3.0 dB ripple) and linear phase filter having the same order.
"""
lpf_flat = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Butterworth',f_start=0.1,f_stop = 4)

lpf_chebyshev = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Chebyshev',f_start=0.1,f_stop = 4,ripple = 3.0)


fig,ax = plt.subplots(2,1)
ax_2 = ax[0].twinx()

lpf_flat.plot_summary()
#lpf_chebyshev.plot_summary(ax = ((ax[0],ax_2),ax[1]))

#ax[0].set_ylim(-20,0)
#ax_2.set_ylim(-40,0)


class test():
    def __init__(self,a,b):
        self.a = a
        self.b = b
    
    def summation(self):
        return self.a + self.b
