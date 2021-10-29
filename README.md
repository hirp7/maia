# maia, microwave passive filter design assistant 
is a support tool for RF passive filter design based on microwave filter synthesis theory.
Using the insertion loss method, one can obtain a required filter design parameters such as g-parameters, transformed inductance and capacitance for ladder LC circuit from the princple low pass filter in accordance with filter requirements such as center frequency fc, fractional bandwidth FBW, ripple band Lar, and the number of stage N, 
A class *LowPassFilter()* is easily created with filter specification as shown below.

    lpf_flat = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Butterworth',f_start=0.1,f_stop = 4)
    lpf_chebyshev = LowPassFilter(N=5,fc=2.0e9,Base='C_base',Type='Chebyshev',f_start=0.1,f_stop = 4,ripple = 3.0)

The latest version of maia covers following filter-types:

- Lowpass Filter (LPF)
- Bandpass Filter (BPF)

in either butterworth or chebyshev characteristic.

Then you can either see the cartesian graph of filter characteristics, S-parameter and group delay by a method *.plot_summary()*

![lpf_example_S](/images/lpf_example_S.png)


or the schematic of LC ladder circuits with corresponding inductance and capacitance by a method *.schematic()*
![lpf_example_S](/images/lpf_example_schematic.png)




