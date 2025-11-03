# Companion repository to "Empirical plunge profiles of time-frequency localization operators"

<div align="center"> <img width="90%" src="figures/blob_a10_M100.png"/> </div>

This repository contains code to repeat the main experiments of the paper <a href="https://arxiv.org/abs/2502.11805">Empirical plunge profiles of time-frequency localization operators</a>. Inside the `LTFAT` folder, <a href="https://github.com/SimonHalvdansson/Time-Frequency-Plunge-Profiles/blob/main/ltfat-2.6.0/plunge_profile.m">`plunge_profile.m`</a> contains the code for the main figures in the paper. Change the `load_symbol` argument to change the symbol. The code for the RGB figure can be found in <a href="https://github.com/SimonHalvdansson/Time-Frequency-Plunge-Profiles/blob/main/ltfat-2.6.0/acc_spec_color.m">`acc_spec_color.m`</a> file. For the figure comparing different windows, use the `acc_spec_color.m` routine.

The python file <a href="https://github.com/SimonHalvdansson/Time-Frequency-Plunge-Profiles/blob/main/f(k)_plot.py">`f(k)_plot.py`</a> contains a numerical verification of Proposition 2.2 which was used to produce Figure 1 in the paper.

### Abstract
For time-frequency localization operators, related to the short-time Fourier transform, with symbol $R\Omega$, we work out the exact large $R$ eigenvalue behavior for rotationally invariant $\Omega$ and conjecture that the same relation holds for all scaled symbols $R \Omega$ as long as the window is the standard Gaussian. Specifically, we conjecture that the $k$-th eigenvalue of the localization operator with symbol $R\Omega$ converges to $\frac{1}{2}erfc\big( \sqrt{2\pi}\frac{k-R^2|\Omega|}{R|\partial \Omega|} \big)$ as $R \to \infty$. To support the conjecture, we compute the eigenvalues of discrete frame multipliers with various symbols using LTFAT and find that they agree with the behavior of the conjecture to a large degree.</div>

### Sorted eigenvalues validation
<div align="center"> <img width="90%" src="figures/sorted.png"/> </div>

### Window dependence on eigenvalue spectrogram separation 
<div align="center"> <img width="50%" src="figures/rgb_specs.png"/> </div>
