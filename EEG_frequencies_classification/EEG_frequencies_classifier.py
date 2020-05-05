import pyedflib
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss

class EEG_frequencies_classifier:
    
    def __init__(self):
        self.t = None #time
        self.s = None #signal
        self.signal_labels = None
        self.Fs = None #frequency sampling
        self.freq = None
        self.Pxx = None
        self.chosen_channels = None
        self.time_range = None
    
    def read_from_edf(self, filename, Fs): 
        print('Reading the data from edf file.')
        f = pyedflib.EdfReader(filename)
        self.Fs = Fs
        n = f.signals_in_file
        signal_labels = f.getSignalLabels()
        self.s = np.zeros((n, f.getNSamples()[0]))
        self.signal_labels = {}
        for i in np.arange(n):
            self.s[i, :] = f.readSignal(i)
            self.signal_labels[signal_labels[i]] = i 
        self.t = np.arange(0, len(self.s[0,:])/Fs, 1/Fs)
        
        
    def montage(self):
        '''common_averge'''
        print('Making common average montage of the signal.')
        self.s -= np.mean(self.s,0)
        
        
    def filter_signal(self, n, Wn, btype):
        ''' Filter the signal with Butterworth filter.
        n : int
        The order of the filter.
        Wn : float or length-2 numpy array
        The critical frequency or frequencies in Hz. For lowpass and highpass filters 
        Wn is a scalar; for bandpass and bandstop filters Wn is a length-2 numpy array.
        btype : {‘lowpass’, ‘highpass’, ‘bandpass’, ‘bandstop’}
        The type of the filter.
        '''
        print('Filtering the signal with {} Hz {} Butterworth filter of order {}.'.format(Wn, btype, n))
        if type(Wn) == list:
            Wn = np.array(Wn)
        b, a = ss.butter(n, Wn/(self.Fs/2), btype = btype)
        self.s = ss.filtfilt(b, a, self.s)
               
                
    def channels_to_str(self, channels):
        '''if channels is a list of channels' numbers
        converts channels' numbers to channels' names
        else return unchanged channels'''
        if type(channels[0]) == int:
            return [list(self.signal_labels.keys())[i] for i in channels]
        return channels
   
    def channels_to_num(self, channels):
        '''if channels is a list of channels' names
        converts channels' names to channels' numbers
        else return unchanged channels'''
        if type(channels[0]) == str:
             return [self.signal_labels[i] for i in channels]
        return channels
    
    def power_spectrum(self, channels, time_range=None):
        '''channels : list of channels numbers or names;
        time_range : length-2 list of the begin and the end of signal fragment in seconds,
        default is full time of signal duration
        '''
        print('Computing power spectrum.')
        if time_range == None:
            time_range=[self.t[0],self.t[-1]]
        self.time_range = time_range
        time_range = np.array(time_range)*self.Fs
        time_range = time_range.astype(int)
        self.chosen_channels = self.channels_to_num(channels)
        signal_fragment = self.s[self.chosen_channels, time_range[0]:time_range[1]]
        self.freq, self.Pxx = ss.welch(signal_fragment, self.Fs)
        return self.freq, self.Pxx
        
    def plot_signal(self, channels, time_range=None):
        '''channels : list of channels numbers or names;
        time : length-2 list of the begin and the end of signal fragment in seconds,
        default is full time of signal duration
        '''
        if time_range == None:
            time_range=[self.t[0],self.t[-1]]
        self.time_range = time_range
        time_range = np.array(time_range)*self.Fs
        time_range = time_range.astype(int)
        signal_fragment = self.s[self.channels_to_num(channels), time_range[0]:time_range[1]]
        chosen_channels = self.channels_to_str(channels)
        
        plt.figure(figsize=(9,6.5))
        for i, ch in enumerate(chosen_channels):
            plt.subplot(len(chosen_channels), 1, i+1)
            plt.subplots_adjust(hspace=0.1)
            plt.suptitle('Fragment of EEG signal')
            plt.plot(self.t[time_range[0]:time_range[1]], signal_fragment[i,:])
            plt.ylabel(ch+'\n $[µV]$', fontsize=11) #plt.ylabel(ch, fontsize=12)

            ax = plt.gca()
            ax.yaxis.set_label_coords(-0.05, 0.5)
            plt.yticks(fontsize=9, rotation=0)
            if i != len(chosen_channels)-1:
                plt.xticks([])
            if i == len(chosen_channels)-1:
                plt.xlabel('[s]')
        plt.show()
        
        
    def plot_power_spectrum(self, frequencies=None):
        '''frequencies : length-2 list of the begin and the end of range of frequencies in Hz,
        default is maximum frequencies range for power spectrum;
        plot power spectrum of EEG signal computed in power_spectrum
        for chosen range of frequencies'''
        
        if frequencies == None:
            frequencies = [self.freq[0], self.freq[-1]]
        
        chosen_channels = self.channels_to_str(self.chosen_channels)
        
        plt.figure(figsize=(9,6.5))
        for i, ch in enumerate(chosen_channels):
            plt.subplot(len(chosen_channels), 1, i+1)
            plt.subplots_adjust(hspace=0.1)
            plt.suptitle('Power spectrum of EEG fragment from {} to {} s of signal'.format(self.time_range[0],self.time_range[1]))
            plt.plot(self.freq, self.Pxx[i])
            plt.ylabel(ch+'\n $[µV^{2}]$', fontsize=11) #plt.ylabel(ch, fontsize=12)

            ax = plt.gca()
            ax.yaxis.set_label_coords(-0.05, 0.5)
            plt.yticks(fontsize=9, rotation=0)
            plt.xlim(frequencies)
            if i != len(chosen_channels)-1:
                plt.xticks([])
            if i == len(chosen_channels)-1:
                plt.xlabel('[Hz]')
        plt.show()

        