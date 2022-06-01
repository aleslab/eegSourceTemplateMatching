# TemplateSourceLoc
functions for EEG source localization using templates 

1. Make sure that the templates and your data share the same EEG montage and reference 
(i.e. they have the same electrode locations and the same electrode or average reference). 
If it is not the case, run createCustomTemplates to create templates that fit your data.

2. Use sourceLoc.m to get source activity in ROIs 
This function uses the best regularisation value based on L-curve
In case you get a warning or the regularisation does not seem right: 
a. Try plotting the L-curve when computing sourceLoc: sourceLoc(template , data, 1)
b. Visually determine the corner of the curve.
c. Use the value of this corner in sourceLocFix as a fixed regularisation value: sourceLocFix(template , data, CornerValue)

For more details see help in function
