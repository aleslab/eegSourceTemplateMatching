% Example script 
% step 1. simulates activity over time in V1 in condition 1 and hMT+ in
% condition 2 for a montage with 32 electrodes. [go directly to next steps
% for examples in using the toolbox]
% step 2. the templates corresponding to the 32 electrodes montage are   
% created using createCustomTemplates. (For EGI, please use the EGI
% montages as they cover electrode locations outside the standard 10-05
% system)
% step 3. The activity is recovered from functional 
% brain areas using fitEEGTemplates. 
% step 4. plots the recovered activity in different brain areas
% If it works, you should see some figures showing the alignment of the
% electrodes with the specified 32 electrodes montage. 
% The last 2 figures should look similar. One is the simulated activity,
% the other one is the retrieved activity over time (with noise) in the 18 
% functional brain areas for the 2 conditions (in blue and red).  


%% step 1. simulate activity for 2 conditions 
% (corresponding to the EEG average across participants in 2 conditions)
load('template_Standard_1005.mat')
% use a montage with 32 electrodes
setOfChannels = {'Fp1' 'AF3' 'F7' 'F3' 'FC1' 'FC5' 'T7' 'C3' 'CP1' 'CP5'...
    'P7' 'P3' 'Pz' 'PO3' 'O1' 'Oz' 'O2' 'PO4' 'P4' 'P8' 'CP6' 'CP2' 'C4'...
    'T8' 'FC6' 'FC2' 'F4' 'F8' 'AF4' 'Fp2' 'Fz' 'Cz'};
% get indexes corresponding to this set of channels from the templates
matchIndex = cell2mat(arrayfun(@(x) cellfind(templates.label,setOfChannels{x}),1:length(setOfChannels),'uni',false));

cond1ROI = [1 2]; % bilateral V1 active  
cond2ROI = [17 18]; % bilateral hMT+ active 
% activation over time using half cosine
% 90 timepoints, peak activation at 45  
cosFilt = cos(-pi:pi/45:pi-pi/45);
cosFilt(cosFilt<0) = 0;
simulCond = zeros(size(templates.listROIs,2),length(cosFilt),2);
simulCond(cond1ROI,:,1) = repmat(cosFilt,length(cond1ROI),1);
simulCond(cond2ROI,:,2) = repmat(cosFilt,length(cond1ROI),1);

% add noise (needs different regularisation when no noise)
simulCond = simulCond + randn(size(simulCond))*0.1;

% simulate data over electrodes
dataCond1 = templates.weights(matchIndex,:) * simulCond(:,:,1);
dataCond2 = templates.weights(matchIndex,:) * simulCond(:,:,2);

% reference data to average 
dataCond1 = bsxfun(@minus,dataCond1, mean(dataCond1));
dataCond2 = bsxfun(@minus,dataCond2, mean(dataCond2));


% read channel locations for the simulated activity
% in EEGLAB channel locations can be loaded from a file using readlocs or 
% can look up locations for a set of labelled electrodes. It is recommended 
% to choose MNI to load the locations of your electrodes. 
% EEG.chanlocs = struct('labels', setOfChannels);
% EEG.nbchan = length(setOfChannels); EEG.chaninfo.nosedir = '+X';
% EEG = pop_chanedit( EEG );

% in Fieldtrip can use ft_read_sens
% The locations of the channels can also be created "manually". Below, it
% uses the cartesian coordinates provided by Biosemi (uses spherical head
% model & RAS coordinate system). The format follows fieldtrip format. 
channelInfo.label = setOfChannels;
% order should follow the order of the channel labels
channelInfo.chanpos = [-27	83	-3
-36	76	24
-71	51	-3
-48	59	44
-33	33	74
-78	30	27
-87	0	-3
-63	0	61
-33	-33	74
-78	-30	27
-71	-51	-3
-48	-59	44
0	-63	61
-36	-76	24
-27	-83	-3
0	-87	-3
27	-83	-3
36	-76	24
48	-59	44
71	-51	-3
78	-30	27
33	-33	74
63	0	61
87	0	-3
78	30	27
33	33	74
48	59	44
71	51	-3
36	76	24
27	83	-3
0	63	61
0	0	88];

%% step 2. create templates corresponding to the montage
% the following uses an average reference (0), which match the reference of
% the data
% Try myTemplates = createCustomTemplates(channelInfo) for other options or
% read the help in the function
myTemplates = createCustomTemplates(channelInfo,0); 
% for EEGLAB: myTemplates = createCustomTemplates(EEG.chanlocs,0); 

%% step 3. fit data
% data simulated with 90 timepoints, peak activation half way (45)
% in 2 conditions (one with V1 activation, the other with hMT+)
% put the 2 conditions in the save variable (assumes that the simulated
% data is an average of multiple participants)
averageData(:,:,1) = dataCond1;
averageData(:,:,2) = dataCond2;
% for a cell it would be: averageData{1} = dataCond1; averageData{2} = dataCond2;
areaActive = fitEEGTemplates(averageData,myTemplates);

%% step 4. plot
figure; hold on
loc = [1:9;10:18]; loc = loc(:);
for area=1:length(templates.listROIs)
    subplot(2,9,loc(area)); hold on
    plot(simulCond(area,:,1))
    plot(simulCond(area,:,2))
    title(templates.listROIs(area))
    ylim([-1 1])
    xlabel('time')
end
legend('cond 1 simulates V1', 'cond 2 simulates hMT+','location','best')

figure; hold on
loc = [1:9;10:18]; loc = loc(:);
for area=1:length(templates.listROIs)
    subplot(2,9,loc(area)); hold on
    plot(areaActive(area,:,1))
    plot(areaActive(area,:,2))
    title(myTemplates.listROIs(area))
    ylim([-1 1])
    xlabel('time')
end
legend('cond 1 retrieves V1', 'cond 2 retrieves hMT+','location','best')
