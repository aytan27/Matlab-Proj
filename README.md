# Matlab-Exp
Requires matlab to run the code. 

The stimulus generation program for the dichotic sample discrimination task is dichotic1b.  This is a listening task where different tones are sent to both ears. 
For example in one condition, Right-Quiet, tones 2, 3,4 are presented in the Right ear at 50 dB while tones 1,3,5,7 are presented at 70 dB in the Left ear 

The data analysis program is Datadis_ins which tells you the minimum # of blocks to generate listener JF's (i.e., sample data from listener JF is provided here) lowest RMS value and the decision weighting strategy of the listener. 
Datadis_ins includes a function call efficiency to calculate efficiency measures for weighting and internal noise (loss of effiency due to non-optimal weighting strategy).


There's also a program ewaifofl_ayt1. This code runs simulations of an ideal performer in a narrowband three-tone task where they have to judge if the signal or standard sounds higher in pitch. This model simulation runs 25K trials and is based off of Feth's Envelope-Weighted Average of Instantaneous Frequency pitch model to predict listening performance in frequency discrimination task, with the modification to account for off-frequency listening. This program gives the mean of the signal and standard as well as a figure of the predicted asymmetric pitch weights to give comparison between observed pitch data and the model. 

There's also sample data in the AT_Pitchdata.zip from subject AT who participated in the pitch discrimination experiment using threetone stimuli. Data analysis file D2data11.m outputs 2 figures detailing average threshold, performance and weights that shows it matches the model weights predicted by the EWAIF off-frequency listening model. 

