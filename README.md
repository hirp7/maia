# maia, microwave passive filter design assistant 
is a support tool for RF passive filter design based on the insertion loss method.
Using integrated modules, one can obtain filter design parameters such as g-parameters for LC-ladder circuit, transformed inductance and capacitance LC from the princple LPF, and *network object in scikit-rf*.
Depending on a filter requirement such as center frequency *fc*, fractional bandwidth *FBW* (for BPF), ripple band *Lar* in dB (for chebyshev-type), and the number of stage *N*, a object *LowPassFilter* or *BandPassFilter* is easily defined as shown below:

    lpf_flat = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Butterworth',f_start=0.1,f_stop = 4)

    lpf_chebyshev = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Chebyshev',f_start=0.1,f_stop = 4,ripple = 3.0)

The latest version of maia covers following filters:

- Lowpass Filter (LPF)
- Bandpass Filter (BPF)

in either *butterworth* or *chebyshev* type.

Then you can either see the cartesian graph of filter characteristics, S-parameter and group delay by a method *.plot_summary()*

![lpf_example_S](/images/LPF_Butterworth_2_4.png)


or the schematic of LC ladder circuits with corresponding inductance and capacitance by a method *.schematic()*

![lpf_example_S](/images/LPF_Butterworth_2_4_schematic.png)




