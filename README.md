# EEG Source Template Matching

A toolbox for fitting ERP data to group average topographic templates created from functional areas mapped by fMRI and the topographies calculated using individually defined BEM forward models. 

This toolbox can be used to determine the contribution of functional brain 
sources to the EEG scalp  signal by fitting these topographic templates (which represent 
the scalp activity of 18 visual areas) to the EEG data.
The EEG data should not contained large artefacts (eye-blinks, muscle, etc.) which could make the localisation fails. 

This toolbox is designed to work with both EEGLAB and FieldTrip. This toolbox contains a rudimentary EEGLAB plugin that can be run from the GUI. For this, a "STUDY" design is required. However, for more control use it as a script (see ExampleScript_EEGLAB_study.m)

Overview of steps to use this toolbox:

1. Find the EEG templates that correspond to your EEG montage. The templates and your 
data should share the same electrodes and reference. 
A set of standard EGI templates are provided in the templates folder. 
You can also create, for any EEG montage but EGI, your own templates based on your 
set-up by running createCustomTemplates. This function will extract templates
from a 10-05 system (high density electrode montage) that matches your set of electrodes. 
mytemplates = createCustomTemplates(mychannellocations)
It is recommended to plot mytemplates as a sanity check (set as default for 
EEGLAB data, otherwise plot as you would for any topography). 
You can save the templates and only do this step once for a given EEG montage.

2. Run fitEEGTemplates to obtain the normalised contributions of functional 
brain sources to your EEG data. This function uses a L-curve regularisation. 
In case you get a warning or the regularisation does not seem right: 
a. Try plotting the L-curve when computing fitEEGTemplates: 
fitEEGTemplates(data, templates, 1)
b. Visually determine the corner of the L-curve.
c. Use the value of this corner: 
fitEEGTemplates(data, templates, CornerValue)



Getting Started options:


ExampleScript.m
This file creates a toy 32 channel montage and simulates activity over time in V1 in condition 1 and hMT+ in condition 2.
Custom templates are created for the 32 channel montage. The activity is recovered in the 18 visual areas and plotted for the 2 conditions.


ExampleScript_EEGLAB_study.m
This file contains an example code to use an existing EEGLAB study dataset.  It shows  how to load a study, create a custom template, fit data to the template and plot the results.
You can use this dataset from EEGLAB for testing:  https://sccn.ucsd.edu/eeglab/download/animal_study.zip



EEGLAB Plugin:

To use it as a EEGLAB plugin, the toolbox should be located in the EEGLAB/plugins folder before starting EEGLAB. 
Load a study. If not already computed compute appropriate ERP's for your study for -all- electrodes (important because template uses all electrods). 
Then under STUDY menu choose "create custom template".  This will save a custom template file.  Look at the figures to ensure the template is aligned with your electrode montage (to within about a cm). 
Then select "Fit Visual Area Templates" this should have your template file already chosen.  This should create a plot with the results of the template fit.





For more details see help in functions.
You can report any bug or suggest improvements at:
https://github.com/aleslab/eegSourceTemplateMatching/issues
