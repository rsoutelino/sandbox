from scipy.signal import butter, lfilter


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import freqz

    # Sample rate and desired cutoff frequencies.
    fs      = (1.) / (3600)        # as in a hourly current meter record
    lowcut  = (1.) / (10*24*3600)     # diurnal signal  
    highcut = (1.) / (3*24*3600)     # semi-diurnal signal

    # Plot the frequency response for a few different orders, to illustrate 
    #     the filter design.
    plt.figure(1)
    b, a = butter_bandpass(lowcut, highcut, fs, order=3)
    w, h = freqz(b, a)
    freq = (fs * 0.5 / np.pi) * w
    peri = (1./(freq)) / 3600
    plt.plot(peri, abs(h))
    plt.fill([12, 24, 24, 12, 12], [-0.5, -0.5, 1.5, 1.5, 0.5], '0.7')
    plt.axis([0, 48*2, -0.5, 1.5])
    plt.xlabel('Period [h]')
    plt.ylabel('Gain')

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)

    # Filter a noisy signal.
    T = 30*24*3600 # 1 month
    nsamples = np.round(T * fs)
    t = np.linspace(0, T, nsamples, endpoint=False)

    # hypothetical M2
    aM2 = 2
    fM2 = 2*np.pi / (12.4 * 3600)
    yM2 = aM2 * np.cos( fM2 * t )

    # hypothetical K1
    aK1 = 0.5
    fK1 = 2*np.pi / (23.93 * 3600)
    yK1 = aK1 * np.cos( fK1 * t )

    # hypothetical residual
    ar = 0.4
    fr = 2*np.pi / (4 * 24 * 3600)   # four days
    yr = ar * np.cos( fr * t ) + 0.6 # with a simple offset 

    # hypothetical noise
    yn = np.random.randn(t.size) / 2

    # computing the time-series
    y = yM2 + yK1 + yr + yn

    plt.figure(2)
    plt.clf()
    plt.plot(t/86400, y, 'r', label='Noisy signal', alpha=0.3)

    yf = butter_bandpass_filter(y - y.mean(), lowcut, highcut, fs, order=3)
    yf += y.mean()
    plt.plot(t/86400, yf, 'r', label='Filtered signal', linewidth=2)
    plt.xlabel('Time [days]')
    plt.hlines( [ (aM2+aK1)*-1+y.mean(), (aM2+aK1)+y.mean() ], 0, T/86400, linestyles='--' )
    plt.grid(True)
    plt.legend(loc='upper left')

    plt.show()