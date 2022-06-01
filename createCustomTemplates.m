function customTemplates = createCustomTemplates(channelInfo,varargin)
% function creating ROI templates corresponding to the EEG montage with 
% which the dataset was recorded
% Plot the topographies of the customTemplates for each ROI to check that
% they are similar to what is in the paper. Done automatically for EEGlab
% input data. If it does not seem right, check the coordinate system of
% your data (see below for details)
% usage EEGlab: customTemplates = createCustomTemplates(EEG.chanlocs)
% usage Fieldtrip: customTemplates = createCustomTemplates(cfg.elec)

% input channelInfo: variable containing electrode definition
% For EEGlab users = EEG.chanlocs (obtained from readlocs or using the GUI)
% For better results use the MNI channel location file for electrode 
% position (default since EEGlab 2021)
% https://eeglab.org/tutorials/04_Import/Channel_Locations.html 
% For Fieldtrip users, structure (usually cfg.elec obtained from 
% ft_read_sens) containing
%   elec.elecpos = Nx3 matrix with carthesian (x,y,z) coordinates of each
%                  electrode
%   elec.label   = cell-array of length N with the label of each electrode
%   elec.chanpos = Nx3 matrix with coordinates of each sensor

% output is a structure containing:
%   customTemplates.data = matrix of EEG response for N electrodes x 18 ROI-templates 
%   customTemplates.ROInames = cell-array of the 18 ROI names (for L-left and R-right separately)
%   customTemplates.chanLabels = cell-array of length N with the label of
%   each electrode (same as the dataset)
%   customTemplates.ref = reference used for the templates (0=average)
%   customTemplates.matchedLabels = best matching labels from the templates

% Optional input (has to be in the following order, [] to skip one):
% 1st: Reference of the EEG montage. Default is read from channelInfo for 
% EEGlab. If not found, ask user or can be entered here. The reference 
% should be the channel number or 0 for average reference. 
% usage: customTemplates = createCustomTemplates(EEG.chanlocs,0)
% 2nd: plot the templates 1 (default) or not 0. Only works for EEG.chanlocs 
% input
% 3rd: Coordinate system: 'ALS' or 'RAS'. By default the coordinate system  
% is ALS for EEGlab users, RAS otherwise but can be specified by the user. 
% (see below for more details. Other coordinate system can potentially be 
% added in the program) 
% usage: customTemplates = createCustomTemplates(EEG.chanlocs,[],[],'RAS')

% To do the fitting, need 3-D cartesian coordinates of the electrodes
% corresponding to the EEG montage that is used AND at least 8 channel 
% labels to match. These can be extracted from the dataset or entered 
% manually when prompted if no label is found.
% Also need the reference of the montage

% The ROI templates were created based on 3D electrode coordinates from the 
% fieldtrip standard_1005.elc file constructed by Robert Oostenveld.
% The file can be read with ft_read_sens. 
% The electrode positions are represented in mm in the MNI coordinate 
% system: RAS for first dimension orients towards Right, the 
% second dimension orients towards Anterior, the third dimension orients 
% towards Superior. 
% see also: https://www.fieldtriptoolbox.org/template/electrode/
% for coordinate system: https://www.fieldtriptoolbox.org/faq/coordsys/
% EEGlab uses a different coordinate system: ALS. x is towards the nose, 
% y is towards the left ear, and z towards the vertex. So before finding 
% the best alignment between the dataset and templates, the coordinates 
% need to be flipped accordingly. 

% template assumed to be in the same directory, if it cannot be found, user 
% is prompted to pick the directory manually

% could potentially use alignFiducials.m to align the coordinates
% automatically

if length(channelInfo) == 1 % not eeglab
    eeglabUser = 0; coordsys = 'RAS';chanLabelsData = channelInfo.label;
    plotMap = 0;refIndex=[];
else
    eeglabUser = 1; coordsys = 'ALS';chanLabelsData = {channelInfo.labels};
    plotMap = 1;
    if isempty(channelInfo(1).ref)
        refIndex = [];
    elseif strcmp(channelInfo(1).ref,'average')
        refIndex = 0;
    else
        refIndex = channelInfo(1).ref;
    end
end
if ~isempty(varargin) 
    refIndex = varargin{1}; % use reference entered by user
    if length(varargin) == 2
        plotMap = varargin{2};
    end
    if length(varargin) == 3
        coordsys = varargin{3}; % use coordinate system entered by user
    end
end
fprintf('Assumes %s coordinate system.\n',coordsys)

% try to match electrodes between dataset and template using standard 10-20: Fp1 Fp2 Fz F7 F3 C3 T7 P3 
% P7 Pz O1 Oz O2 P4 P8 T8 C4 F4 F8 Cz
listLabels = {'Fp1' 'Fp2' 'Fz' 'F7' 'F3' 'C3' 'T7' 'P3' 'P7' 'Pz' 'O1' 'Oz' ...
    'O2' 'P4' 'P8' 'T8' 'C4' 'F4' 'F8' 'Cz' };

% get chan indexes that match between the listLabels and the data
matchIndex = arrayfun(@(x) cellfind(chanLabelsData,listLabels{x}),1:length(listLabels),'uni',false);
matchIndex = cell2mat(matchIndex(~cellfun('isempty', matchIndex))); % keep only non-empty indexes

% if data does not contain channel labels
% ask user to enter a set of channel labels with their corresponding
% indexes from standard 10-20
if length(matchIndex)<8
    while length(matchIndex) < 8
        prompt = ['Could not find matching channels. Enter at least 8 channel ' ...
            'numbers from the list below in the SAME order, separated by space. '...
            'If a channel is not in your EEG montage write 0. '...
            'Fp1 Fp2 Fz F7 F3 C3 T7 P3 P7 Pz O1 Oz O2 P4 P8 T8 C4 F4 F8 Cz'];
        dlgtitle = 'Create custom templates';
        answer = inputdlg(prompt,dlgtitle);
        matchIndex = str2num(answer{1});
        if length(matchIndex)~=length(listLabels)
            msgError = errordlg('Please enter as many numbers as the number of channels in the list','Input error','modal');
            uiwait(msgError);
        else
            if sum(~ismember(matchIndex,0:length(chanLabelsData)))>0 % include 0 for non-existing channel
                msgError = errordlg('Please enter valid channel number','Input error','modal');
                uiwait(msgError);
            else
                matchIndex = matchIndex(matchIndex>0);
            end
        end
    end
end
% list of the labels to match
listLabelMatch = chanLabelsData(matchIndex);

% reference of the EEG montage
if isempty(refIndex) 
    prompt = {'Please enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create custom templates';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end
while length(refIndex)~=1 || ~ismember(refIndex,0:length(chanLabelsData)) % include 0 for average reference
    msgError = errordlg('Please enter a valid number for the EEG reference','Input error','modal');
    uiwait(msgError);
    prompt = {'Enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create custom templates';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end
    
% load templates (including elecDef)
templateFile = ['templates' filesep 'template_Standard_1005.mat'];
templateDir = 1;
while ~exist(templateFile) && templateDir~=0
%     uiwait(msgbox('Could not find template_Standard_1005.mat file','Information','modal'));
    templateDir = uigetdir('','Please select the directory containing template_Standard_1005.mat');
    templateFile = [templateDir filesep 'template_Standard_1005.mat'];
end
load(templateFile)

% get the indexes that match with the data (& from the list of channels)
matchIndexROI = cell2mat(arrayfun(@(x) cellfind(elecDef.label,listLabelMatch{x}),1:length(listLabelMatch),'uni',false));

% get electrodes coordinate to match to the templates
% change orientation if necessary
if eeglabUser
    if strcmp(coordsys,'ALS')
        coordToMatch = [-[channelInfo.Y]; [channelInfo.X]; [channelInfo.Z]]';
    else % RAS
        coordToMatch = [[channelInfo.X]; [channelInfo.Y]; [channelInfo.Z]]';
    end
else % fieldtrip
    if strcmp(coordsys,'ALS')
        coordToMatch = [-channelInfo.chanpos(:,2) channelInfo.chanpos(:,1) channelInfo.chanpos(:,3)];
    else % RAS
        coordToMatch = channelInfo.chanpos;
    end
end


% align the montages using the matching electrodes
% Add dimension for "homogeneous coordinates"  for affine fitting.
coordToMatch(:,4) = 1; 
elecDef.chanpos(:,4) = 1; 
% Affine:   templateLoc = customLoc * AffineMtx
% affineMtx = templateLoc/customLoc
affineMtx = elecDef.chanpos(matchIndexROI,:)' / coordToMatch(matchIndex,:)' ;
% Apply transform to all the electrodes of the montage
elecTrans = (affineMtx*coordToMatch')';


%%% Scatter plot 
% check that the new selected electrodes are close to the ones from the template
cmap = hsv(length(matchIndex)); 
figure;scatter3(elecTrans(matchIndex,1),elecTrans(matchIndex,2),elecTrans(matchIndex,3),[],cmap,'filled')
hold on;scatter3(elecDef.chanpos(matchIndexROI,1),elecDef.chanpos(matchIndexROI,2),elecDef.chanpos(matchIndexROI,3),[],cmap)
for nn=1:length(matchIndex)
    text(elecTrans(matchIndex(nn),1),elecTrans(matchIndex(nn),2),elecTrans(matchIndex(nn),3),listLabelMatch{nn})
end
title('Checking alignment')


%%%%%%% Best Match
% Use knnsearch to find the best match between electrodes of the current 
% montage and the template
% returns electrode indexes (same length as nb of chan) & distances
[bestElec, euclideanDist] = knnsearch(elecDef.chanpos(:,1:3),elecTrans(:,1:3));

% create new custom template using the best matching electrodes
matchedTemplate = avMap(bestElec,:); 

figure; scatter3(elecDef.chanpos(bestElec,1),elecDef.chanpos(bestElec,2),elecDef.chanpos(bestElec,3),'filled')
hold on; scatter3(elecTrans(:,1),elecTrans(:,2),elecTrans(:,3))
title('Best matching electrodes')



%%%%%% interpolate the electrodes (exponential weighted)
% pairwise distance between each electrode locations of the current 
% montage and the template   
distances = pdist2(elecDef.chanpos(:,1:3),elecTrans(:,1:3),  'squaredeuclidean');
% weight distances exponentially
weightedDist = exp(-distances/(20^2));
% normalize
normDist = weightedDist ./ repmat(sum(weightedDist),length(weightedDist),1);
% apply interpolation
interpTemplate = (avMap' * normDist)'; 

    

% re-reference the templates
if refIndex == 0
    % average reference: substract average across channels
    rerefTemplates = bsxfun(@minus,matchedTemplate, mean(matchedTemplate)); 
    rerefTemplatesInterp = bsxfun(@minus,interpTemplate, mean(interpTemplate)); 
else
    % substract the reference channel
    rerefTemplates = bsxfun(@minus,matchedTemplate, matchedTemplate(refIndex,:)); 
    rerefTemplatesInterp = bsxfun(@minus,interpTemplate, interpTemplate(refIndex,:)); 
end

% output structure
customTemplates.data = rerefTemplates;
customTemplates.dataInterp = rerefTemplatesInterp;
customTemplates.chanLabels = chanLabelsData';
customTemplates.ROInames = listROIs';
customTemplates.ref = refIndex;
customTemplates.matchedLabels = chanLabels(bestElec);
customTemplates.matchDistances = euclideanDist; % in mm

% plot optional topographies 
if eeglabUser && plotMap == 1
    loc = [1:9;10:18]; loc = loc(:);
    mm = round(max(max(abs(customTemplates.data))),-1);
    figure('position', [200, 200, 1000, 200])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(customTemplates.data(:,roi),channelInfo,'colormap','jet');
        caxis([-mm mm]);title(listROIs(roi));
    end
end  

    