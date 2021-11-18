# maia, microwave analytically integrated-design assistant 
is a support tool for RF circuit design, especially for microwave passive filter design at this moment.

The latest version of maia covers following filters:

- Lowpass Filter (LPF)
- Bandpass Filter (BPF)

in either *butterworth* or *chebyshev* type based on insertion loss method.

Using integrated modules, one can obtain filter design parameters such as g-parameters for LC-ladder circuit, transformed inductance and capacitance LC from the princple LPF, and *network* object in ```scikit-rf```.
Depending on a filter requirement such as the cut-off frequency for LPF or the center of frequency for BPF *fc*, fractional bandwidth *FBW* (for BPF), ripple band *Lar* in dB (for chebyshev-type), and the number of stage *N*, a object *LowPassFilter* or *BandPassFilter* is easily defined as shown below:
```
    lpf_flat = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Butterworth',f_start=0.1,f_stop = 4)
    lpf_chebyshev = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Chebyshev',f_start=0.1,f_stop = 4,ripple = 3.0)
```

Once a object is defined, one can see the summary of results (Fig. 1) and schematic (Fig. 2) by methods, *.plot_summary()* and *.schematic()*, respectively.


![lpf_example_S](/images/LPF_Butterworth_2_4.png)

(Fig.1)

![lpf_example_S](/images/LPF_Butterworth_2_4_schematic.png)

(Fig.2)



