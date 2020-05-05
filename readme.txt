Project description

EEG_frequencies_classification

Program for computing and plotting power spectrum of EEG signal

Author: Zuzanna Lewandowska, University of Warsaw

EEG_frequencies_classification is a program that enables computing 
the power spectrum of EEG signal for chosen EEG channels in chosen 
range of time and plot it for chosen range of frequencies. EEG signal
should be in EDF format.

Program has been written in Python 3.6.9.

Data used to write a program comes from the PhysioNet databases:
Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, 
Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, 
and PhysioNet: Components of a New Research Resource for Complex 
Physiologic Signals (2003). Circulation. 101(23):e215-e220.
(https://physionet.org/content/eegmmidb/1.0.0/, access: 2020-05-05)
from experiment:
Schalk, G., McFarland, D.J., Hinterberger, T., Birbaumer, N., Wolpaw, 
J.R. BCI2000: A General-Purpose Brain-Computer Interface (BCI) System. 
IEEE Transactions on Biomedical Engineering 51(6):1034-1043, 2004.
(https://www.ncbi.nlm.nih.gov/pubmed/15188875, access: 2020-05-05)

Usage example

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