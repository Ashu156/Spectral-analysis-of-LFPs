# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 12:15:27 2021

@author: Dell
"""
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

# Import filter function
from neurodsp.filt import filter_signal

# Import utilities for loading and plotting data
# from neurodsp.utils.download import load_ndsp_data
from neurodsp.plts.time_series import plot_time_series

##################### Load data ############################
rawdata = scipy.io.loadmat(file_path)
data = rawdata['lfp_14']    # data from a single channel
data = (data).T;
fs = 1000;                  # sampling frequency in Hz
times = rawdata['tx_14']    # time vector
times = (times).T

print(f"DATA: (n_trials, n_times)={data.shape}; SAMPLING FREQUENCY={fs}Hz; "
      f"TIME VECTOR: n_times={len(times)}")

# Define a frequency range to filter the data
sig = np.mean(data, axis  = 0)
f_range = (2, 12)

# Bandpass filter the data, across the band of interest
sig_filt = filter_signal(sig, fs, 'bandpass', f_range, filter_type='iir', butterworth_order = 4)

# Plot filtered signal
# fig = plt.figure(figsize=(14, 6))
plot_time_series(times, [sig, sig_filt], ['Raw', 'Filtered'])
# plt.show()



fig = plt.figure(figsize=(14, 6))
# adding the mean PSD over trials
# plt.subplot(1, 2, 1)
plt.plot(times, sig, color = 'black')
plt.plot(times, sig_filt, color = 'red', linewidth = 5.0)
plt.ylim([-50, 50])
# plt.xlim([-1000, 0])
plt.show()
fig.savefig('Control dmPFC short theta.eps', dpi = 300)

