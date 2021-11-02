MAIA, Microwave pAssive fIlter design Assistant, is a support tool for filter design based on the theory of microwave filter synthesis.
Depending on filter specification such as center frequency, fractional bandwidth, the number of stages, ripple loss, as input parameters,
MAIA accordingly generates g-parameters, lumped (LC) circuit, and S-parameter Network file.


MAIA currently covers following filter-types:
-Lowpass Filter (LPF)
-Bandpass Filter (BPF)
in either butterworth or chebyshev characteristic.


Using a method .plot_summary() S-parameter results as well as group delay characteristics are summarized in a figure.
Generated S-parameter is directly related to class Network in scikit-rf, so that they are integrated in further RF Analysis.
A Schematic view is also generated using Circuit-drawing library,SchematicDraw().