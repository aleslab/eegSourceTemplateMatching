% Example script 
% This script simulates activity over time in V1 in condition 1 and hMT+ in
% condition 2 for a montage with 16 electrodes. 
% Then, templates corresponding to the 16 electrodes montage are created  
% using createCustomTemplates. The activity is recovered from functional 
% brain areas using fitEEGTemplates and plotted for the 2 conditions
% separately. 
% If it works, you should see some figures showing the alignment of the
% electrodes with the specified 16 electrodes montage and the retrieved
% activity in the 18 functional brain areas (over time). Note that some
% activity is retrieved in non-active areas probably due to the very small
% number of electrodes used here (16) but the max activity is retrieved in
% the correct areas. 

% for cellfind to work outside the function (or when running the script one
% line at a time), needs to add it to path! 

%% simulate activity for 2 conditions 
% (corresponding to the EEG average across participants in 2 conditions)
load(['templates' filesep 'template_Standard_1005.mat'])
% use a montage with 16 electrodes
setOfChannels = {'Fp1' 'Fp2' 'F4' 'Fz' 'F3' 'T7' 'C3' 'Cz'...
    'C4' 'T8' 'P4' 'Pz' 'P3' 'O1' 'Oz' 'O2'};
% get indexes corresponding to this set of channels from the templates
matchIndex = cell2mat(arrayfun(@(x) cellfind(templates1005.label,setOfChannels{x}),1:length(setOfChannels),'uni',false));

% topography for activation in bilateral V1  
topo1 = sum(templates1005.weights(matchIndex,[1 2]),2); 
% topography for activation in bilateral hMT+
topo2 = sum(templates1005.weights(matchIndex,[17 18]),2); 
% activation over time using half cosine
% 90 timepoints, peak activation at 45  
cosFilt = cos(-pi:pi/45:pi-pi/45);
cosFilt(cosFilt<0) = 0;
dataCond1 = topo1 .* cosFilt;
dataCond2 = topo2 .* cosFilt;

% get simulated data over time for plotting
simulatedData = zeros(size(templates1005.listROIs,2),length(cosFilt),2);
simulatedData(1,:,1) = cosFilt;
simulatedData(2,:,1) = cosFilt;
simulatedData(17,:,2) = cosFilt;
simulatedData(18,:,2) = cosFilt;

% read channel locations for the simulated activity
% in EEGLAB can be loaded from a file using readlocs or can look up 
% locations for a set of labelled electrodes. It is recommended to choose 
% MNI to load the locations of your electrodes. Here it will be:
% chanlocs = struct('labels', setOfChannels);
% channelInfo = pop_chanedit( chanlocs );

% in Fieldtrip can use ft_read_sens
% The locations of the channels can also be created "manually". Below, it
% uses the cartesian coordinates provided by Biosemi (uses spherical head
% model & RAS coordinate system). The format follows fieldtrip format. 
channelInfo.label = setOfChannels;
% order should follow the order of the channel labels
channelInfo.chanpos = [-27	83	-3
    27	83	-3
    48	59	44
    0	63	61
    -48	59	44
    -87	0	-3
    -63	0	61
    0	0	88
    63	0	61
    87	0	-3
    48	-59	44
    0	-63	61
    -48	-59	44
    -27	-83	-3
    0	-87	-3
    27	-83	-3];

%% create templates corresponding to the montage
% the following uses an average reference. 
% Try myTemplates = createCustomTemplates(channelInfo) for other options or
% read the help in the function
myTemplates = createCustomTemplates(channelInfo,0); 

%% fit data
% data simulated with 90 timepoints, peak activation half way (45)
% in 2 conditions (one with V1 activation, the other with hMT+)
% put the 2 conditions in the save variable (simulation is the ERP of an
% average of multiple participants)
averageData(:,:,1) = dataCond1;
averageData(:,:,2) = dataCond2;
% for a cell it would be: averageData{1} = dataCond1; averageData{2} = dataCond2;
areaActive = fitEEGTemplates(averageData,myTemplates);

%% plot
figure; hold on
loc = [1:9;10:18]; loc = loc(:);
for area=1:size(areaActive,1)
    subplot(2,9,loc(area)); hold on
    plot(simulatedData(area,:,1))
    plot(simulatedData(area,:,2))
    title(myTemplates.listROIs(area))
    ylim([-1 1])
    xlabel('time')
end
legend('cond 1 simulates V1', 'cond 2 simulates hMT+','location','best')

figure; hold on
loc = [1:9;10:18]; loc = loc(:);
for area=1:size(areaActive,1)
    subplot(2,9,loc(area)); hold on
    plot(areaActive(area,:,1))
    plot(areaActive(area,:,2))
    title(myTemplates.listROIs(area))
    ylim([-1 1])
    xlabel('time')
end
legend('cond 1 retrives V1', 'cond 2 retrieves hMT+','location','best')
