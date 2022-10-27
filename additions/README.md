# The 18 ROIs were created from the MRI and fMRI scans of 50 individuals.
There is also the possibility to create 20 ROIs which separate V1 into V1
ventral and V1 dorsal. However, these 2 ROIs are based on 27 individuals. 
The templates for these 20 ROIs are in this folder. 
Just use them directly for EGI system. 
To create your own, use 
myTemplates = createCustomTemplates(channelInfo,[],[],[],[],1)