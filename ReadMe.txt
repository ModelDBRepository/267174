
This is the readme for the models associated with the paper:

Mauro Ursino*, Giulia Ricci, Laura Astolfi, Floriana Pichiorri, Manuela Petti and Elisa Magosso, (2021)
A Novel Method to Assess Motor Cortex Connectivity and Event Related Desynchronization Based on Mass Models
Brain sciences


The patient data used in the study are not published in the repository as the patient has not given explicit consent in line with the current GDPR.
For this reason, to run the codes it is first necessary to create a .mat file, named “Spettri_dati.mat”, that contains the Power Spectral Densities 
(PSD) and the Coherences (C) of the experimental signals. 

The original patient EEG data have been acquired with Fc=100 Hz in three different experimental conditions: baseline, affected and unaffected. 
Then, cortical sources have been reconstructed and 6 cortical regions of interest (ROI) have been identified. 

Data in “Spettri_dati.mat” must be organised as follows:
number of frequency points (nSample): 501
number of cortical Regions of Interest (nROI): 6

P_baseline_media, P_affected_media, P_unaffected_media:  are the matrices of Power Spectral Densities of six cortical regions of interest (nROI) 
in three different experimental conditions. Specifically, each matrix size must be: [nSample x nROI]  

C_baseline_media, C_affected_media, C_unaffected_media: are the {nROI x nROI} cell arrays of coherences (Cij) in the three experimental conditions 
(baseline, affected and unaffected). Specifically, each cell array element contains the coherence between two cortical regions (ROIi and ROIj) 
and has dimension: [nSample x 1]. Each row (i) of the cell is characterized by 6 columns of elements (j=1, 2, ...6), each containing a specific 
coherence between the i-th ROI and the j-th ROI (for example, the cell elements on the first row (i=1) of the cell array are: C11, C12, C13, C14, C15, C16). 

Once all spectra and coherences have been computed for the 3 experimental conditions, it is necessary to save all matrices and cell arrays in a file 
named “Spettri_dati.mat”. The file must then be copied into both sub-folders (programmi34e56 and programmi1e2) of the main folder (programmi base), 
as at the beginning of the codes the file “Spettri_dati.mat” is loaded. 





This main folder (programmi base) contains two sub-folders.


FIRST FOLDER NAME: programmi34e56
Contains the programs used for parameter estimation in SMAp L, SMAp R, PMD L and PMD R regions, and for subsequent simulations.

“Stima_4ROI.m”: the program loads at line 7 the experimental data, contained in “Spettri_dati.mat”.
 Then, it executes the parameter estimation starting from the parameter values contained in the file at line 105 
(in the present example the file is named “prova_4ROI_2giugno2020.mat”). The program makes use of the function “costo_4ROI.m” (at line 123).
 At the end of the program, by uncommenting line 600, the estimated parameters can be saved (the present file name “prova_4ROI” can be modified”). 

“costo_4ROI.m”: this function computes the cost function, starting from the parameter values given as argument. 
 To this end, it in turn calls a second function named “modello_fitting.m”.

“modello_fitting.m”: this function integrates the set of model differential equations. 
It receives as input the parameters and provides as output the quantities of the model (state variables and activity of neuronal populations).

The other functions “outfun”, “outfung” or “outfunp” are used for the minimization procedure. 
At present, the program uses “patternsearch” and the function “outfunp” (defined at line 144) in the options of this method, 
but other methods can be used by uncommenting or commenting lines 134-155.

“Simulazione_4ROI”: the program performs the simulation and plots all figures (Power Spectra Densities and Coherences). At first, the program loads the experimental data, contained in “Spettri_dati.mat” (loaded at line 6). Then, parameter values are loaded at line 56 (in the present example these are contained in the file “prova_4ROI_2giugno2020.mat”). Moreover, by uncommenting lines 474 and 475, the simulated activity of neural populations and the EEG signals can be saved (examples are given in the file “uscite_4ROI_2giugno2020.mat” and “eeg_4ROI_2giugno2020.mat”).


 

SECOND FOLDER NAME: programmi1e2
Contains the programs used for parameters estimations in the regions M1h L and M1h R, and for subsequent simulations.

“Stima_12_inib.m”: the program loads the experimental data at line 8, contained in “Spettri_dati.mat”.
Moreover, it loads the signals feed-forwarded to the regions M1hL and M1hR from the upstream four ROIs (SMAp L, SMAp R, PMD L, PMD R) 
at line 9 (in the present example the file name is “uscite_4ROI_2giugno2020.mat”). 
Then, it executes the parameter estimation starting from the parameter values contained in the file loaded at line 77 
(in the present example the file is named “prova12inib_19giugno2020.mat”). The program makes use of the function “costo1e2_inib.m” (line 98).
 At the end of the program, by uncommenting the line 345, the parameters can be saved in a file (the present name “prova” can be changed”). 

“costo1e2_inib.m”: this function computes the cost function, starting from the parameter values given as argument.  

The other functions “outfun”, “outfung” or “outfunp” are used for the minimization procedure. 
At present, the program uses “patternsearch” and the function “outfunp” (defined at line 111) in the options of this method,
 but other methods can be used by uncommenting or commenting lines 102-122.

“Simulazione1e2_inib”: the program performs the simulation and plots all figures (Power Spectra Densities and Coherences).
 At first, it loads the experimental data, contained in “Spettri_dati.mat” (at line 5). 
Moreover, it loads the simulated activities of neural populations and the EEG signals feed-forwarded from the upstream four ROIs 
(SMApL, SMApR, PMD L, PMD R) previously computed with the program “Simulazione_4ROI.m”. 
In the present example these files are named “uscite_4ROI_2giugno2020.mat” (loaded at line 6) and “eeg_4ROI_2giugno2020.mat” (loaded at line 7). 
The parameter values are loaded at line 90 (in the present example these are contained in the file “prova12inib_19giugno2020.mat”). 
By uncommenting lines 643 and 644, the simulated activities in the populations and the EEG signals can be saved 
(examples are given in the file “uscite_12_inib.mat” and “eeg_12_inib.mat”). 

Finally, the same files but without the suffix “inib” (i.e., “Stima_12.m”, “costo1e2.m” and “Simulazione1e2.m”) can be used in exactly 
the same way as the other programs, but assuming the presence of an excitatory connection between each SMAp and the contralateral M1h 
(see Supplementary Material III). 



