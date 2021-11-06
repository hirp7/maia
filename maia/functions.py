
import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,brentq,minimize, differential_evolution
from skrf.media import DefinedGammaZ0,MLine,Freespace
import schemdraw
import schemdraw.elements as elm
import copy


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


from schemdraw.segments import *
        
class TransmissionLine(elm.Element2Term):
    def __init__(self, *d, **kwargs):
        super().__init__(*d, **kwargs)
        radius = 0.075
        width = 3
        length = 3
        self.segments.append(Segment([(0, 0), (length, 0)],lw = 4.0))
        self.segments.append(Segment([(0, -width), (length, -width)],lw = 4.0))

        self.segments.append(SegmentCircle((0, 0), radius,fill = None))
        self.segments.append(SegmentCircle((0, -width), radius,fill = None))
        self.segments.append(SegmentCircle((length, -width), radius,fill = None))
        self.segments.append(SegmentCircle((length, 0), radius,fill = None))

        self.endpoints((0,0),(length,0))


        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (0, -width)
        self.anchors['p3'] = (length, 0)
        self.anchors['p4'] = (length,-width)

class Two_Port_Box(elm.Element):
    def __init__(self,length, *d, **kwargs):
        super().__init__(*d, **kwargs)
        BoxWidth = 4
        BoxLength = 2
        width = 3
        #length = 0.5
        self.segments.append(Segment([(0, -width), (length, -width)],lw = 1))
        
        #Box
        self.segments.append(Segment([(length, -width/2 - BoxWidth/2), (length+BoxLength, -width/2 - BoxWidth/2)],lw = 2))
        self.segments.append(Segment([(length, -width/2 - BoxWidth/2), (length, -width/2 + BoxWidth/2)],lw = 2))
        self.segments.append(Segment([(length+BoxLength, -width/2 - BoxWidth/2), (length+BoxLength, -width/2 + BoxWidth/2)],lw = 2))
        self.segments.append(Segment([(length, -width/2 + BoxWidth/2), (length+BoxLength, -width/2 + BoxWidth/2)],lw = 2))


        self.segments.append(Segment([(length + BoxLength, 0), (length+BoxLength+length, 0)],lw = 1))
        self.segments.append(Segment([(length + BoxLength, -width), (length+BoxLength+length, -width)],lw = 1))
        self.segments.append(Segment([(0, 0), (length, 0)],lw = 1))

        
        """
        self.segments.append(Segment([(0, -width), (length, -width)],lw = 4.0))
        self.segments.append(SegmentCircle((0, 0), radius,fill = None))
        self.segments.append(SegmentCircle((length, 0), radius,fill = None))
        self.segments.append(SegmentCircle((0, -width), radius,fill = None))
        self.segments.append(SegmentCircle((length, -width), radius,fill = None))
        """
        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (0, -width)
        self.anchors['p3'] = (length*2+BoxLength, 0)
        self.anchors['p4'] = (length*2+BoxLength,-width)
        self.anchors['p1'] = (0, 0)
      



g_butterworth = {
"1":[2.0000,1.0000],
"2":[1.4142,1.4142,1.0000],
"3":[1.0000,2.0000,1.0000,1.0000],
"4":[0.7654,1.8478,1.8478, 0.7654, 1.0000],
"5":[0.6180, 1.6180, 2.0000, 1.6180, 0.6180, 1.0000],
"6":[0.5176, 1.4142, 1.9318, 1.9318, 1.4142, 0.5176, 1.0000],
"7":[0.4450, 1.2470, 1.8019, 2.0000, 1.8019, 1.2470, 0.4450,1.0000],
"8":[0.3902, 1.1111, 1.6629, 1.9615, 1.9615, 1.6629, 1.1111, 0.3902, 1.0000],
"9":[0.3473, 1.0000, 1.5321, 1.8794, 2.0000, 1.8794, 1.5321, 1.0000, 0.3473, 1.0000],
"10":[0.3129, 0.9080, 1.4142, 1.7820, 1.9754, 1.9754, 1.7820, 1.4142, 0.9080, 0.3129, 1.0000]
}




"""
g = {
     'Butterworth':g_butterworth,
     'Chebyshev':g_chebyshev1,
     'Chebyshev':g_chebyshev,
     'PhaseLinear':g_phaselinear}
"""

def g_butterworth(N):
    g = [1]
    def g_n(n):
        return np.round(2*sin((2*n-1)*pi/2/N),4)
    g_add = list(map(g_n,np.arange(N)+1))
    g.extend(g_add)
    g.append(1.0)
    return g


def g_chebyshev(N,ripple):
    g = [1]
    #ripple = 10**(ripple/10)
    beta = np.log(1/np.tanh(ripple/17.37))
    gamma = np.sinh(beta/2/N)
    def g_n(n):
        if n==1:
            return 2/gamma*sin(pi/2/N)
        else:
            #return 2*sin((2*g_n(n-1)*pi/2/N))
            return 1/g_n(n-1)*4*sin((2*n-1)*pi/2/N)*sin((2*n-3)*pi/2/N)/\
                (gamma**2 + sin((n-1)*pi/N)**2)
            
    g_add = list(map(g_n,np.arange(N)+1))
    g.extend(g_add)
    if N%2==1:
        g.append(1.0)
    else:
        g.append(1/np.tanh(beta/4)**2)
    return np.round(g,4)

def matrix_conversion(L): #LC list is converted to numpy matrix mxnxn
    n = len(L[0])
    return np.array([L[0],L[1]]).reshape(n,2)

def matrix_conversion2(L): #LC list is converted to numpy matrix mxnxn
    n = len(L)
    return np.array(L).reshape(n,1)

#mxnxn
def unit_change(A):
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
    


def LPF(g,fc=1e9,filter_type = 'L_base',Freq = rf.Frequency(1,10,1001,'ghz')): #Freq is defined rf.Frequency
    global LC
    Omega_c = 2*pi*fc
    media_FS = Freespace(frequency = Freq)
    L_Value = []  #attain L and C value
    C_Value = []
    #load = g[-1]
    #g = g[:-1]
        
    if filter_type == 'L_base':
        list(map(lambda x:L_Value.append(50*x/Omega_c),g[::2]))
        list(map(lambda x:C_Value.append(x/50/Omega_c),g[1::2]))
        LC = [None]*len(L_Value+C_Value)
        LC[::2] = list(map(media_FS.inductor,L_Value))
        LC[1::2] = list(map(media_FS.shunt_capacitor,C_Value))
            
    elif filter_type == 'C_base':
        list(map(lambda x:L_Value.append(50*x/Omega_c),g[1::2]))
        list(map(lambda x:C_Value.append(x/50/Omega_c),g[::2]))
        LC = [None]*len(L_Value+C_Value)
        LC[1::2] = list(map(media_FS.inductor,L_Value))
        LC[::2] = list(map(media_FS.shunt_capacitor,C_Value))
    
    def N_cascade(Networks): #Networks is a list containing N networks
        length = len(Networks)
        cascaded = Networks[0]
        for i in range(length-1):
            cascaded = cascaded**Networks[i+1]
        return cascaded

    Filter = N_cascade(LC)
    Filter.renormalize(50)
    L = matrix_conversion2(L_Value)
    C = matrix_conversion2(C_Value)

    return Filter,[L,C]#np.array([L_Value,C_Value]),LC  

def BPF(g,fc,FBW,filter_type = 'L_base',Freq = rf.Frequency(1,10,1001,'ghz')):
    #f1 = fc - fc*FBW/2
    #f2 = fc + fc*FBW/2
    Omega_c = 2*pi*fc
    #FBW = 2*pi*(f2 - f1)/Omega_c
    media_FS = Freespace(frequency = Freq)
    #load = g[-1]
    #g = g[:-1]
    L_Value = []  
    L_C_Value = []
    C_Value = []
    C_L_Value = []
    LC_BPF = []
    
    if filter_type == 'L_base':
        list(map(lambda x:L_Value.append(50*x/Omega_c/FBW),g[::2]))
        list(map(lambda x:L_C_Value.append(1/x/50/Omega_c*FBW),g[::2]))
        
        list(map(lambda x:C_Value.append(x/50/Omega_c/FBW),g[1::2]))
        list(map(lambda x:C_L_Value.append(1/x*50/Omega_c*FBW),g[1::2]))

        LC1 = [None]*len(L_Value+C_Value)
        LC2 = [None]*len(L_Value+C_Value)
        LC1[::2] = list(map(media_FS.inductor,L_Value))
        LC1[1::2] = list(map(media_FS.shunt_capacitor,C_Value))
        LC2[::2] = list(map(media_FS.capacitor,L_C_Value))
        LC2[1::2] = list(map(media_FS.shunt_inductor,C_L_Value))

            
    elif filter_type == 'C_base':
        list(map(lambda x:L_Value.append(50*x/Omega_c/FBW),g[1::2]))
        list(map(lambda x:L_C_Value.append(1/x/50/Omega_c*FBW),g[1::2]))
        
        list(map(lambda x:C_Value.append(x/50/Omega_c/FBW),g[::2]))
        list(map(lambda x:C_L_Value.append(1/x*50/Omega_c*FBW),g[::2]))

        LC1 = [None]*len(L_Value+C_Value)
        LC2 = [None]*len(L_Value+C_Value)
        LC1[1::2] = list(map(media_FS.inductor,L_Value))
        LC1[::2] = list(map(media_FS.capacitor,C_Value))
        LC2[1::2] = list(map(media_FS.shunt_capacitor,L_C_Value))
        LC2[::2] = list(map(media_FS.shunt_inductor,C_L_Value))
                
    def N_cascade(Networks1,Networks2): #Networks is a list containing N networks
        length = len(Networks1)
        cascaded = Networks1[0]**Networks2[0]
        for i in range(length-1):
            cascaded = cascaded**Networks1[i+1]**Networks2[i+1]
        return cascaded
                
    Filter = N_cascade(LC1,LC2)
    Filter.renormalize(50)  
    L = matrix_conversion([L_Value,L_C_Value])
    C = matrix_conversion([C_Value,C_L_Value])
    return Filter,[L,C]

def Schematic_LPF(LC,base ='L_base'):
    d = schemdraw.Drawing() #L-Base or C-base odd or even    
    if base == 'L_base':
        n = len(LC[1])
        Start = d.add(elm.Resistor().label('$R_S$ = 50 Ohm'))
                
        def LC_Schematic(L,C,d = d):
            d += elm.Inductor().label(str(np.round(unit_change(L)[0][0],2)) + unit_change(L)[1] +'H') 
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(C)[0][0],2)) + unit_change(C)[1] +'F')  # Go off in another direction temporarily
            #print('d.here:', d.here)
            d.pop()   # Return to the pushed position/direction
        
        list(map(LC_Schematic,LC[0],LC[1]))
        
        if len(LC[0]) > n:
            d += elm.Inductor().label(str(np.round(unit_change(LC[0][-1])[0][0],2)) + unit_change(LC[0][-1])[1] +'H')
        else:
            d += elm.Line()
        
    if base == 'C_base':
        n = len(LC[0])

        Start = d.add(elm.Resistor().label('$R_S$ = 50 Ohm'))
                
        def CL_Schematic(L,C,d = d): #LC is a list
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(C)[0][0],2)) + unit_change(C)[1] +'F')  # Go off in another direction temporarily
            #print('d.here:', d.here)
            d.pop()   # Return to the pushed position/direction
            d += elm.Inductor().label(str(np.round(unit_change(L)[0][0],2)) + unit_change(L)[1] +'H')

        
        list(map(CL_Schematic,LC[0],LC[1]))
        
        if len(LC[1]) > n:                       
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(LC[1][-1])[0][0],2)) + unit_change(LC[1][-1])[1] +'F')             
            d.pop()   # Return to the pushed position/direction
            d += elm.Line()
        
    d.add(elm.Resistor().down().label('$R_L$ = 50 Ohm',rotate = True,loc = 'bottom'))    
    d.add(elm.Line().left().tox(Start.start))
        
    d.add(elm.SourceSin().up())
    d.draw()


        #if it´s odd there is not
    # Now we want to close the loop, but can use `tox`
    # to avoid having to know exactly how far to go.
    # Note we passed the [x, y] position of capacitor C,
    # but only the x value will be used.
    return d



#LC = [[[1,2,3],[2,3,4]],[[3,4,5],[4,5,6]]]

def Schematic_BPF(LC,base = 'L_base'):
    d = schemdraw.Drawing()    
    
    if base == 'L_base':
        n = len(LC[1]) #the number of capacitor
        Start = d.add(elm.Resistor().label('$R_S$ = 50 Ohm'))        
        def LC_Schematic(L,C,d = d):
            d += elm.Inductor().label(str(np.round(unit_change(L[0])[0],2)) + unit_change(L[0])[1] +'H')
            d += elm.Capacitor().label(str(np.round(unit_change(L[1])[0],2)) + unit_change(L[1])[1] +'F')
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(C[0])[0],2)) + unit_change(C[0])[1] +'F')  # Go off in another direction temporarily
            d.pop()   # Return to the pushed position/direction
            d += elm.Line()
            d.push()
            d += elm.Inductor().down().label(str(np.round(unit_change(C[1])[0],2)) + unit_change(C[1])[1] +'H')  # Go off in another direction temporarily            
            d.pop()
        list(map(LC_Schematic,LC[0],LC[1]))
        #a = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]] b = [x[1] for x in a] list comprehension
        if len(LC[0]) > n:
            d += elm.Inductor().label(str(np.round(unit_change(LC[0][0][-1])[0],2)) + unit_change(LC[0][0][-1])[1] +'H')
            d += elm.Capacitor().label(str(np.round(unit_change(LC[0][1][-1])[0],2)) + unit_change(LC[0][1][-1])[1][-1] +'F')
        else:
            d += elm.Line()

    if base == 'C_base':
        n = len(LC[0]) #number of inductor
        Start = d.add(elm.Resistor().label('$R_S$ = 50 Ohm'))        
        def CL_Schematic(L,C,d = d):
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(C[0])[0],2)) + unit_change(C[0])[1] +'F')  # Go off in another direction temporarily
            d.pop()   # Return to the pushed position/direction
            d += elm.Line()
            d.push()
            d += elm.Inductor().down().label(str(np.round(unit_change(C[1])[0],2)) + unit_change(C[1])[1] +'H')  # Go off in another direction temporarily            
            d.pop()

            d += elm.Inductor().label(str(np.round(unit_change(L[0])[0],2)) + unit_change(L[0])[1] +'H')
            d += elm.Capacitor().label(str(np.round(unit_change(L[1])[0],2)) + unit_change(L[1])[1] +'F')

        list(map(CL_Schematic,LC[0],LC[1]))

        if len(LC[1]) > n:
            d.push()  # Save this drawing position/direction for later
            d += elm.Capacitor().down().label(str(np.round(unit_change(LC[1][0][-1])[0],2)) + unit_change(LC[1][0][-1])[1] +'F')  # Go off in another direction temporarily
            d.pop()   # Return to the pushed position/direction
            d += elm.Line()
            d.push()
            d += elm.Inductor().down().label(str(np.round(unit_change(LC[1][1][-1])[0],2)) + unit_change(LC[1][1][-1])[1] +'H')  # Go off in another direction temporarily            
            d.pop()
            d += elm.Line()


            
        
        #d += elm.Capacitor() #if it´s odd there is not
    d.add(elm.Resistor().down().label('$R_L$ = 50 Ohm',rotate = True,loc = 'bottom'))    
    # Now we want to close the loop, but can use `tox`
    # to avoid having to know exactly how far to go.
    # Note we passed the [x, y] position of capacitor C,
    # but only the x value will be used.
    d.add(elm.Line().left().tox(Start.start))
    
    d.add(elm.SourceSin().up())
    #d.draw()
    return d


"""

def impedance(w,Z0,h,t,ep_r,tanD,Freq,fc): #to find a required impedance, Freq is rf.Frequency, fc is float
    return MLine(frequency = Freq[str(fc)],w = w*1e-3, h=h,t = t,ep_r = ep_r, tand=tanD).Z0-Z0


def StepImpedanceLPF(N,g,Z0,h,t,ep_r,tanD,Freq,fc):
    #N and len(g) must be equal
    #Z0 = [Source, High, Low] #Source, High, Low
    #'C-base', 'L-base'
    w_source = np.round(brentq(impedance,0.001,20,args = (Z0[0],h,t,ep_r,tanD,Freq,fc)),2) #w_minimum,w_max,
    w_1 = np.round(brentq(impedance,0.001,20,args = (Z0[1],h,t,ep_r,tanD,Freq,fc)),2) #w_minimum,w_max,
    w_2 = np.round(brentq(impedance,0.001,20,args = (Z0[2],h,t,ep_r,tanD,Freq,fc)),2) #w_minimum,w_max,
    w = [w_source,w_1,w_2]


    beta = []
    
    for i in range(len(w)):
        beta.append(MLine(frequency = Freq[str(fc)],w = w[i]*1e-3,h=h,t = t,ep_r = ep_r,\
                          tand = tanD).beta)
            
    length = []
    
    
    Line = []#Thru-Step-Thru
    Line.append(MLine(frequency = Freq,w = w[0]*1e-3,h=h,t = t,ep_r = ep_r,\
                          tand = tanD).line(5e-3,'m',embed=True,z0 = 50)) #Load


    for i in range(len(g)):
        if i%2==0:
            length.append(g[i]*Z0[2]/50/beta[2]) #high impedance
            Line.append(MLine(frequency = Freq,w = w[2]*1e-3,h=h,t = t,ep_r = ep_r,\
                      tand = tanD).line(length[i],'m'))

        else:
            length.append(g[i]*50/Z0[1]/beta[1]) #low impedance
            Line.append(MLine(frequency = Freq,w = w[1]*1e-3,h=h,t = t,ep_r = ep_r,\
                      tand = tanD).line(length[i],'m'))
        
    
    Line.append(MLine(frequency = Freq,w = w[0]*1e-3,h=h,t = t,ep_r = ep_r,\
              tand = tanD).line(5e-3,'m',embed=True,z0 = 50)) #Source

    
    S = Line[0]
    for i in range(len(Line)-1):
        S = S ** Line[i+1]
        
    return S,length,Line,w,beta

"""
            
def Step_Schematic(Z,length):
    d = schemdraw.Drawing()
    C = d.add(elm.Resistor().label(r'$Z_{source}$'))
    for i in range(len(length)):
        if i%2==0:
            d += TransmissionLine().label(r'$Z_{0(%d)} = $'%(i+1)+str(Z[2])+r' [$\Omega]$',fontsize=12)\
                .label(r'$l_%d = $'%(i+1) + str(np.round(length[i][0]*1e3,2)),loc='bottom',fontsize=12)
        else:
            d += TransmissionLine().label(r'$Z_{0(%d)} = $'%(i+1)+str(Z[1])+r' [$\Omega]$',fontsize=12)\
                .label(r'$l_%d = $'%(i+1) + str(np.round(length[i][0]*1e3,2)),loc='bottom',fontsize=12)
            
    d += elm.Resistor().down().label(r'$Z_{load}$',rotate = True,loc = 'bottom')
    d += elm.Line().left().tox(C.start)
    d.draw()
    return d





