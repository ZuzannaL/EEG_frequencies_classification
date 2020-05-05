from EEG_frequencies_classifier import EEG_frequencies_classifier

filename = r'../data/S001R01.edf' 
EEG = EEG_frequencies_classifier() #init an object
EEG.read_from_edf(filename, Fs=160) #load data from the edf file, frequency sampling of data is 160
EEG.montage() #making montage of the signal
EEG.filter_signal(n = 2, Wn = [49.5, 50.5], btype = 'bandstop') #filtering 
EEG.filter_signal(n = 2, Wn = 1, btype = 'highpass')
print(max(EEG.t)) #checking the duration of signal
EEG.plot_signal([8,9,11,50], [20,30]) #plotting for channels chosen by number from 20 to 30 s of signal
print(EEG.signal_labels) #checking channels' labels
EEG.plot_signal(['C4..','Cz..','C3..','O2..','O1..']) #plotting full signal for channels chosen by name
EEG.power_spectrum([1,4,5,10], [20,30]) #computing power spectrum for channels chosen by number from 20 to 30 s of signal
EEG.plot_power_spectrum() #plotting computed power spectrum for full range of freqencies
EEG.plot_power_spectrum([8,18]) #plotting computed power spectrum for frequencies from 8 to 18 Hz