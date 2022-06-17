% Example script 

%% simulate activity for 2 conditions 
% (corresponding to the EEG average across participants in 2 conditions)
load(['templates' filesep 'template_Standard_1005.mat'])
% use a montage with 16 electrodes
setOfChannels = {'Fp1' 'Fp2' 'F4' 'Fz' 'F3' 'T7' 'C3' 'Cz'...
    'C4' 'T8' 'P4' 'Pz' 'P3' 'O1' 'Oz' 'O2'};
% get indexes corresponding to this set of channels from the templates
matchIndex = cell2mat(arrayfun(@(x) cellfind(templates1005.label,setOfChannels{x}),1:length(setOfChannels),'uni',false));

% topography for activation in bilateral MT  
topo1 = sum(templates1005.weights(matchIndex,[17 18]),2); 
% topography for activation in bilateral V1 
topo2 = sum(templates1005.weights(matchIndex,[1 2]),2); 
% activation over time using half cosine
% 90 timepoints = 90, peak activation at 45  
cosFilt = cos(-pi:pi/45:pi-pi/45);
cosFilt(cosFilt<0) = 0;
dataCond1 = topo1 .* cosFilt;
dataCond2 = topo2 .* cosFilt;


%% read channel locations
% in EEGLAB can look up locations for a set of labelled electrodes
% locations can also be loaded from a file using readlocs
% in Fieldtrip can use ft_read_sens 

if exist('pop_chanedit.m','file') == 2
    % here just use 16 electrodes
    chanlocs = struct('labels', setOfChannels);
    channelInfo = pop_chanedit( chanlocs );
    % in the pop-up window, preferably choose: ?Use MNI coordinates for the BEM
    % Dipfit model? which is what has been used to create the templates.
    % However, for visuallisation try BESA.
    
else
    % Locations of the channels can also be created "manually"
    % use the cartesian coordinates provided by Biosemi (uses spherical head &
    % RAS coordinate system). The format follows fieldtrip format
    channelInfo.label = setOfChannels;
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
end

%% create templates corresponding to the montage
myTemplates = createCustomTemplates(channelInfo);

%% fit data
% data simulated with 90 timepoints = 90, peak activation half way (45)
% in 2 conditions (one with V1 activation, the other with hMT+)
% put the 2 conditions in the save variable (simulation is the ERP of an
% average of multiple participants)
averageData(:,:,1) = dataCond1;
averageData(:,:,2) = dataCond2;
% for a cell it would be: averageData{1} = dataCond1; averageData{2} = dataCond2;
areaActive = fitEEGTemplates(averageData,myTemplates);

%% plot
figure; hold on
for area=1:size(areaActive,1)
    subplot(9,2,area); hold on
    plot(areaActive(area,:,1))
    plot(areaActive(area,:,2))
    title(myTemplates.listROIs(area))
    ylim([-1 1])
end
