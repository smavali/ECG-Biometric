from scipy.signal import butter, filtfilt, iirnotch
import numpy as np

def butter_lowpass(cutoff, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_highpass(cutoff, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_bandpass(lowcut, highcut, sample_rate, order=2):
    nyq = 0.5 * sample_rate
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
    
def filter_signal(data, cutoff, sample_rate, order=2, filtertype='lowpass', return_top = False):
    if filtertype.lower() == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    elif filtertype.lower() == 'notch':
        b, a = iirnotch(cutoff, Q = 30, fs = sample_rate)
    else:
        raise ValueError('filtertype: %s is unknown, available are: \
lowpass, highpass, bandpass, and notch' %filtertype)

    filtered_data = filtfilt(b, a, data)
    
    if return_top:
        return np.clip(filtered_data, a_min = 0, a_max = None)
    else:
        return filtered_data

def remove_baseline_wander(data, sample_rate, cutoff=0.05):
        return filter_signal(data = data, cutoff = cutoff, sample_rate = sample_rate,filtertype='notch')

# >>> from scipy import signal
# >>> import matplotlib.pyplot as plt
# >>> fs = 200.0  # Sample frequency (Hz)
# >>> f0 = 60.0  # Frequency to be removed from signal (Hz)
# >>> Q = 30.0  # Quality factor
# >>> # Design notch filter
# >>> b, a = signal.iirnotch(f0, Q, fs)
# >>> # Frequency response
# >>> freq, h = signal.freqz(b, a, fs=fs)
# >>> # Plot
# >>> fig, ax = plt.subplots(2, 1, figsize=(8, 6))
# >>> ax[0].plot(freq, 20*np.log10(abs(h)), color='blue')
# >>> ax[0].set_title("Frequency Response")
# >>> ax[0].set_ylabel("Amplitude (dB)", color='blue')
# >>> ax[0].set_xlim([0, 100])
# >>> ax[0].set_ylim([-25, 10])
# >>> ax[0].grid()
# >>> ax[1].plot(freq, np.unwrap(np.angle(h))*180/np.pi, color='green')
# >>> ax[1].set_ylabel("Angle (degrees)", color='green')
# >>> ax[1].set_xlabel("Frequency (Hz)")
# >>> ax[1].set_xlim([0, 100])
# >>> ax[1].set_yticks([-90, -60, -30, 0, 30, 60, 90])
# >>> ax[1].set_ylim([-90, 90])
# >>> ax[1].grid()
# >>> plt.show()