# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:14:26 2021

@author: Dell
"""
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()

import numpy as np
import scipy.io

from tensorpac.utils import  PSD

import matplotlib.pyplot as plt

##################### Load data ############################
rawdata = scipy.io.loadmat(file_path)
data = rawdata['lfp_67']    # data from a single channel
data = (data).T;
sf = 1000;                  # sampling frequency in Hz
times = rawdata['tx_67']    # time vector
times = (times).T

print(f"DATA: (n_trials, n_times)={data.shape}; SAMPLING FREQUENCY={sf}Hz; "
      f"TIME VECTOR: n_times={len(times)}")

# function for adding the sections rest / stimulus / post-stimulus to each figure
def add_motor_condition(y_text, fontsize=14, color='k', ax=None):
    x_times = [-1000, 0, 1300]
    x_conditions = ['REST', 'STIMULUS\nON', 'STIMULUS\nOFF']
    if ax is None: ax = plt.gca()  # noqa
    plt.sca(ax)
    plt.axvline(0., lw=2, color=color)
    plt.axvline(1.5, lw=2, color=color)
    for x_t, t_t in zip(x_times, x_conditions):
        plt.text(x_t, y_text, t_t, color=color, fontsize=fontsize, ha='center',
                 va='center', fontweight='bold')

# define time indices where rest, planning and execution are defined
time_rest = slice(500, 1000)
time_prep = slice(1500, 2800)
time_exec = slice(2800, 4000)


data_rest = data[..., time_rest]
psd = PSD(data_rest, sf)

fig = plt.figure(figsize=(14, 6))
# adding the mean PSD over trials
plt.subplot(1, 2, 1)
ax = psd.plot(confidence = None, f_min = 1, f_max=30, log=True, grid=True)
plt.axvline(2, lw=2, color='red')
plt.axvline(12, lw=2, color='red')


# adding the single trial PSD
plt.subplot(1, 2, 2)
psd.plot_st_psd(cmap='gray_r', f_min=1, f_max=80, vmax=3e4, vmin=0., log=False,
                grid=True)
plt.axvline(2, lw=2, color='red')
plt.axvline(12, lw=2, color='red')
plt.tight_layout()
plt.show()
fig.savefig('Single trial PSD baseline.eps', dpi = 300)

data_prep = data[..., time_prep]
psd = PSD(data_prep, sf)

fig = plt.figure(figsize=(14, 6))
# adding the mean PSD over trials
plt.subplot(1, 2, 1)
ax = psd.plot(confidence = None, f_min = 1, f_max=30, log=True, grid=True)
plt.axvline(2, lw=2, color='red')
plt.axvline(12, lw=2, color='red')

# adding the single trial PSD
plt.subplot(1, 2, 2)
psd.plot_st_psd(cmap='gray_r', f_min=1, f_max=80, vmax=3e4, vmin=0., log=False,
                grid=True)
plt.axvline(2, lw=2, color='red')
plt.axvline(12, lw=2, color='red')
plt.tight_layout()
plt.show()
fig.savefig('Single trial PSD stimulus.eps', dpi = 300)