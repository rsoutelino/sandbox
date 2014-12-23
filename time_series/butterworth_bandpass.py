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
    """
    Simple IIR filter
    This is just a use case to run as a script
    Play with lowT and highT parameters to define de band pass filter
    The default parameters consists of a hypothetical time series
        with M2 and K1 tides added to a 4-day residual cosine and 
        some offset elevation, added to a random noise. There are two 
        filter examples, one to preserve M2 and other to preserve k1. 
        Plots illustrate the filter functioning.  

    """

    # creating a hypothetical time series with M2, K1, residual with offset and noise
    plt.figure(figsize=(14,15), facecolor='w')

    fs = 1. / 3600  # as in a hourly current meter record

    T = 30*24*3600 # 1 month
    nsamples = np.round(T * fs)
    t = np.linspace(0, T, nsamples, endpoint=False)

    # hypothetical M2
    aM2 = 1
    fM2 = 2*np.pi / (12.4 * 3600)
    yM2 = aM2 * np.cos( fM2 * t )

    # hypothetical K1
    aK1 = 0.2
    fK1 = 2*np.pi / (23.93 * 3600)
    yK1 = aK1 * np.cos( fK1 * t )

    # hypothetical residual
    ar = 0.2
    fr = 2*np.pi / (4 * 24 * 3600)   # four days
    yr = ar * np.cos( fr * t ) + 0.3 # with a simple offset 

    # hypothetical noise
    yn = np.random.randn(t.size) / 4

    # computing the time-series with the full signal
    y = yM2 + yK1 + yr + yn


    # example for isolating K1 ------------------------------------------------------

    # Sample rate and desired cutoff frequencies.
    lowT, highT = 22, 26               # in hours, to isolate the diurnal signal
    lowcut  = 1. / (highT*3600)         
    highcut = 1. / (lowT*3600)        

    # Plot the frequency response for a few different orders, to illustrate 
    #     the filter design.
    plt.subplot(411)
    b, a = butter_bandpass(lowcut, highcut, fs, order=3)
    w, h = freqz(b, a)
    freq = (fs * 0.5 / np.pi) * w
    peri = (1./(freq)) / 3600
    plt.plot(peri, abs(h), 'k', linewidth=3)
    plt.fill([lowT, highT, highT, lowT, lowT], [-0.5, -0.5, 1.5, 1.5, 0.5], '0.9')
    plt.axis([0, 60, -0.5, 1.5])
    plt.xlabel('Period [h]')
    plt.ylabel('Gain')

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)])
    plt.xlabel('Period [h]')
    plt.ylabel('Gain')
    plt.grid(True)

    # filtering the series and plotting the results
    plt.subplot(412)
    plt.plot(t/86400, y, 'k', label='Noisy signal', alpha=0.3)

    yf = butter_bandpass_filter(y - y.mean(), lowcut, highcut, fs, order=3)
    yf += y.mean()
    plt.plot(t/86400, yf, 'k', label='band-pass for K1', linewidth=2)
    plt.xlabel('Time [days]')
    plt.hlines( [ aK1*-1+y.mean(), aK1+y.mean() ], 0, T/86400, linestyles='--' )
    plt.grid(True)
    plt.legend(loc='upper right', prop={'size': 10})


    # example for isolating M2 ------------------------------------------------------

    # Sample rate and desired cutoff frequencies.
    lowT, highT = 10, 14               # in hours, to isolate the diurnal signal
    lowcut  = 1. / (highT*3600)         
    highcut = 1. / (lowT*3600)        

    # Plot the frequency response for a few different orders, to illustrate 
    #     the filter design.
    plt.subplot(413)
    b, a = butter_bandpass(lowcut, highcut, fs, order=3)
    w, h = freqz(b, a)
    freq = (fs * 0.5 / np.pi) * w
    peri = (1./(freq)) / 3600
    plt.plot(peri, abs(h), 'r', linewidth=3)
    plt.fill([lowT, highT, highT, lowT, lowT], 
             [-0.5, -0.5, 1.5, 1.5, 0.5], 'r', alpha=0.2)
    plt.axis([0, 60, -0.5, 1.5])
    plt.xlabel('Period [h]')
    plt.ylabel('Gain')

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)])
    plt.xlabel('Period [h]')
    plt.ylabel('Gain')
    plt.grid(True)

    # filtering the series and plotting the results
    plt.subplot(414)
    plt.plot(t/86400, y, 'k', label='Noisy signal', alpha=0.3)

    yf = butter_bandpass_filter(y - y.mean(), lowcut, highcut, fs, order=3)
    yf += y.mean()
    plt.plot(t/86400, yf, 'r', label='band-pass for M2', linewidth=2)
    plt.xlabel('Time [days]')
    plt.hlines( [ aM2*-1+y.mean(), aM2+y.mean() ], 0, T/86400, linestyles='--' )
    plt.grid(True)
    plt.legend(loc='upper right', prop={'size': 10})

    plt.show()