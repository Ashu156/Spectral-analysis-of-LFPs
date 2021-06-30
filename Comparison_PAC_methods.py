# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:12:31 2021

@author: Dell
"""
from tensorpac import Pac
from tensorpac.signals import pac_signals_wavelet

import matplotlib.pyplot as plt

f_pha = 10      # frequency phase for the coupling
f_amp = 100     # frequency amplitude for the coupling
n_epochs = 20   # number of trials
n_times = 4000  # number of time points
sf = 512.       # sampling frequency
data, time = pac_signals_wavelet(sf=sf, f_pha=f_pha, f_amp=f_amp, noise=.8,
                              n_epochs=n_epochs, n_times=n_times)

# define a pac object and extract high-resolution phases and amplitudes using
# Morlet's wavelets
p = Pac(f_pha='hres', f_amp='hres', dcomplex='wavelet')

# Extract all of the phases and amplitudes
phases = p.filter(sf, data, ftype='phase', n_jobs=1)
amplitudes = p.filter(sf, data, ftype='amplitude', n_jobs=1)

plt.figure(figsize=(14, 8))
for i, k in enumerate([1, 2, 3, 4, 5, 6]):
    # switch method of PAC
    p.idpac = (k, 0, 0)
    # compute only the pac without filtering
    xpac = p.fit(phases, amplitudes)
    # plot
    plt.subplot(2, 3, k)
    title = p.method.replace(' (', f' ({k})\n(')
    p.comodulogram(xpac.mean(-1), title=title, cmap='viridis')
    if k <= 3:
        plt.xlabel('')

plt.tight_layout()
plt.show()