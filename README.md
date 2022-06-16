# EEG Source Template Matching
Determine the normalised contribution (from -1 to 1) of functional brain 
sources to the EEG scalp  signal by fitting EEG templates (which represent 
the scalp activity of 18 visual areas) to the EEG data.

1. Find the EEG templates that correspond to your EEG montage. The templates and your 
data should share the same electrodes and reference. 
A set of standard EGI templates are provided in the templates folder. 
You can also create, for any EEG montage, your own templates based on your 
set-up by running createCustomTemplates. This function will extract templates
from a 10-05 system (high density electrode montage) that matches your set of electrodes. 
mytemplates = createCustomTemplates(mychannellocations)
It is recommended to plot mytemplates as a sanity check (set as default for 
EEGLAB data, otherwise plot as you would for any topography). 
You only need to do this step once for a given EEG montage.

2. Run fitEEGTemplates to obtain the normalised contributions of functional 
brain sources to your EEG data. This function uses a L-curve regularisation. 
In case you get a warning or the regularisation does not seem right: 
a. Try plotting the L-curve when computing fitEEGTemplates: 
fitEEGTemplates(data, templates, 1)
b. Visually determine the corner of the L-curve.
c. Use the value of this corner: 
fitEEGTemplates(data, templates, CornerValue)

For more details see help in functions.
You can report any bug or suggest improvements at:
https://github.com/aleslab/eegSourceTemplateMatching/issues
