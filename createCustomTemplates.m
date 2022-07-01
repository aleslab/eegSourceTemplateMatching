function customTemplates = createCustomTemplates(channelInfo,varargin)
%createCustomTemplates creates EEG templates according to user's montage. 
%   CUSTOMTEMPLATES = createCustomTemplates(channelInfo) returns
%   topographic EEG templates for 18 functional brain regions according to
%   the channel locations defined in channelInfo. The templates can be saved 
%   and used for other experiments with the same montage.   
%   CUSTOMTEMPLATES is a structure defined with the following fields:
%   	'weights'               - Scalp response. N-by-M matrix with N electrodes and 
%                                 M brain regions (18). 
%       'listROIs'              - Names of the 18 brain regions (for L: left and R: right 
%                                 separately).  
%   	'label'                 - Cell-array of length N with the electrode label from
%                                 the dataset
%       'reference'             - Reference used to create the templates (0=average; otherwise channel number)
%       'matchedLabel'          - Best matching labels from the 10-05 system
%   	'matchedDistances'      - Euclidean distances (in mm) between electrodes
%                                 comparing user & standard EEG montage  
%   	'electrodesIncludedIndex'	- Index of electrodes included in the custom templates
%       'electrodesExcludedIndex'   - Index of electrodes excluded from
%                                     the user's montage (due to missing location information 
%                                     or no satisfactory matching)
%       'electrodesExcludedLabels'	- Labels of excluded electrodes  
%   It is recommended to plot the topographies (1 per brain region) returned from createCustomTemplates (done 
%   automatically for EEGLAB data). If it does not seem right, check the 
%   coordinate system of your data.
%
%   CHANNELINFO should contain 3-D cartesian (x,y,z) coordinates of the electrodes 
%   corresponding to the user's EEG montage. If CHANNELINFO does not contain 
%   at least 8 channel labels that match with a standard 10-20 montage, the
%   user is asked to manually enter a set of matching electrodes. 
%   For EEGLAB users CHANNELINFO is EEG.chanlocs, obtained from readlocs.m or 
%   when adding channel locations from the GUI. For better results use the 
%   MNI channel location file (default since EEGLAB 2021)
%   For Fieldtrip users, CHANNELINFO is cfg.elec obtained from ft_read_sens 
%   which should contain at least elec.label and elec.chanpos
%
%   OPTIONAL INPUT (has to be entered in the following order, use [] to skip one)
%   customTemplates = createCustomTemplates(channelInfo [, reference] [, plot] [, coordsys] [, interpol]) 
%       'reference' - Reference of the EEG montage. If not specified, user 
%                     is prompted to enter it. It should be a channel index (number) 
%                     or 0 for average reference.
%       'plot'      - 1 or 0; 1 plots the templates. Default is 1 for
%                     EEGLAB. This option does not work for fieldtrip
%       'coordsys'  - Coordinate system of the electrode locations: 'ALS' or 'RAS'. 
%                     By default the coordinate system is ALS for EEGLAB 
%                     users, RAS otherwise. Other coordinate system can 
%                     potentially be added in future release.   
%       'interpol'  - default 0, set to 1 to use interpolation method 
%                     instead of closest electrode to fit the montage
%   
%   USAGE: 
%       % Create custom EEG templates from data analysed with EEGLAB
%       mytemplates = createCustomTemplates(EEG.chanlocs)
%       % Create custom EEG templates from data analysed with fieldtrip and use an interpolation method to fit the montages
%       mytemplates = createCustomTemplates(cfg.elec,[],[],[],1)

% Copyright (C) 2022 Marlene Poncet & Justin Ales, University of St
% Andrews, marlene.poncet@gmail.com, jma23@st-andrews.ac.uk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% addpath('subfunctions')

% 4 maximum optional inputs
if nargin>4
    error('Maximum 4 optional inputs')
end

if length(channelInfo) == 1 % not eeglab
    eeglabUser = 0; coordsys = 'RAS';chanLabelsData = channelInfo.label;
    % Check if any electrodes with empty locations
    % Don't think chanpos can ever be empty with fieldtrip struct...
    tf = arrayfun(@(k) ~isempty(channelInfo.chanpos(k)), 1:length(channelInfo.chanpos));
    elecExcludedIndex = find(~tf); elecIncludedIndex = find(tf);
else
    eeglabUser = 1; coordsys = 'ALS';chanLabelsData = {channelInfo.labels};
    % Check if any electrodes with empty locations
    tf = arrayfun(@(k) ~isempty(channelInfo(k).X), 1:numel(channelInfo));
    elecExcludedIndex = find(~tf); elecIncludedIndex = find(tf);
end
if ~isempty(elecExcludedIndex)
    fprintf('%d electrodes do not have location data \n',sum(~tf))
end

% deal with optional inputs
% first set defaults
optargs = {[],1,coordsys,0};
% find empty user inputs
indexOpt = cellfun(@(x) ~isempty(x), varargin);
% overwrite defaults skipping empty ones
optargs(indexOpt) = varargin(indexOpt);
% put in usable variables
[refIndex, plotMap, coordsys, interpolation] = optargs{:};
fprintf('Assumes %s coordinate system.\n',coordsys)

    
% try to match electrodes between dataset and template using standard 10-20: Fp1 Fp2 Fz F7 F3 C3 T7 P3 
% P7 Pz O1 Oz O2 P4 P8 T8 C4 F4 F8 Cz
listLabels = {'Fp1' 'Fp2' 'Fz' 'F7' 'F3' 'C3' 'T7' 'P3' 'P7' 'Pz' 'O1' 'Oz' ...
    'O2' 'P4' 'P8' 'T8' 'C4' 'F4' 'F8' 'Cz' };


    
% get chan indexes that match between the listLabels and the data
matchIndex = arrayfun(@(x) cellfind(chanLabelsData(elecIncludedIndex),listLabels{x}),1:length(listLabels),'uni',false);
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
                listLabelMatch = listLabels(matchIndex>0);
                matchIndex = matchIndex(matchIndex>0);
            end
        end
    end
else
    % list of the labels to match
    listLabelMatch = chanLabelsData(elecIncludedIndex(matchIndex));
end


% reference of the EEG montage
if isempty(refIndex) 
    prompt = {'Please enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create custom templates';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end


% if EEGLAB 'average" label is input set to refIndex to 0; 
if ischar(refIndex) && strmatch(lower(refIndex),'average')
    refIndex = 0;
end

while length(refIndex)~=1 || ~ismember(refIndex,0:length(chanLabelsData)) % include 0 for average reference
    msgError = errordlg('Please enter a valid number for the EEG reference','Input error','modal');
    uiwait(msgError);
    prompt = {'Enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create custom templates';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end
    
% load standard 10-05 templates 
% File assumed to be in the same directory, if it cannot be found, user 
% is prompted to pick the directory manually
templateFile = ['templates' filesep 'template_Standard_1005.mat'];
templateDir = 1;
while ~exist(templateFile,'file') && templateDir~=0
%     uiwait(msgbox('Could not find template_Standard_1005.mat file','Information','modal'));
    templateDir = uigetdir('','Please select the directory containing template_Standard_1005.mat');
    templateFile = [templateDir filesep 'template_Standard_1005.mat'];
end
load(templateFile)

% get the indexes that match with the data (& from the list of channels)
matchIndexROI = cell2mat(arrayfun(@(x) cellfind(templates1005.label,listLabelMatch{x}),1:length(listLabelMatch),'uni',false));


% The templates were created based on 3D electrode coordinates from the 
% fieldtrip standard_1005.elc file constructed by Robert Oostenveld.
% The file can be read with ft_read_sens. 
% The electrode positions are represented in mm in the MNI coordinate 
% system: RAS for first dimension orients towards Right, the 
% second dimension orients towards Anterior, the third dimension orients 
% towards Superior. 
% see also: https://www.fieldtriptoolbox.org/template/electrode/
% and: https://www.fieldtriptoolbox.org/faq/coordsys/
% EEGLAB uses a different coordinate system: ALS. x is towards the nose, 
% y is towards the left ear, and z towards the vertex. 
% So before finding the best alignment between the dataset and templates, 
% the coordinates need to be flipped accordingly.
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
templates1005.chanpos(:,4) = 1; 
% Affine:   templateLoc = customLoc * AffineMtx
% affineMtx = templateLoc/customLoc
affineMtx = templates1005.chanpos(matchIndexROI,:)' / coordToMatch(matchIndex,:)' ;
% Apply transform to all the electrodes of the montage
elecTrans = (affineMtx*coordToMatch')';


%%% Scatter plot 
% check that the new selected electrodes are close to the ones from the template
cmap = hsv(length(matchIndex)); 
figure;scatter3(elecTrans(matchIndex,1),elecTrans(matchIndex,2),elecTrans(matchIndex,3),[],cmap,'filled')
hold on;scatter3(templates1005.chanpos(matchIndexROI,1),templates1005.chanpos(matchIndexROI,2),templates1005.chanpos(matchIndexROI,3),[],cmap)
for nn=1:length(matchIndex)
    text(elecTrans(matchIndex(nn),1),elecTrans(matchIndex(nn),2),elecTrans(matchIndex(nn),3),listLabelMatch{nn})
end
title('Checking alignment')


%%%%%%% Best Match
% Use knnsearch to find the best match between electrodes of the current 
% montage and the template
% returns electrode indexes (same length as nb of chan) & distances
[bestElec, euclideanDist] = knnsearch(templates1005.chanpos(:,1:3),elecTrans(:,1:3));

% detect electrodes that are very far from any existing electrodes in the
% template. These are probably not in the template so should be removed
% from the data (which would be the case for electrodes on the face or
% below the Nz/T9/Iz/T10 line as in EGI montage).
badMatchingElec = find(euclideanDist>=30); % in mm
keepElec = setdiff(1:length(bestElec),badMatchingElec);
if badMatchingElec>0
    disp(['Could not find satisfactory match for electrode(s): ' num2str(elecIncludedIndex(badMatchingElec))])
    disp('One reason could be that these electrodes are not present in the 10-05 standard system.')
    elecExcludedIndex = [elecExcludedIndex elecIncludedIndex(badMatchingElec)];
    elecIncludedIndex = setdiff(elecIncludedIndex, elecIncludedIndex(badMatchingElec));
end

% create new custom template
if ~interpolation
    % create new custom template using the closest electrodes
    % BUT without the bad matching electrodes
    matchedTemplate = templates1005.weights(bestElec(keepElec),:);
else
    % interpolate the electrodes (exponential weighted)
    % interpolation is done including far off electrodes which are exluded 
    % only at the end. 
    % pairwise distance between each electrode locations of the current
    % montage and the template
    distances = pdist2(templates1005.chanpos(:,1:3),elecTrans(:,1:3),  'squaredeuclidean');
    % weight distances exponentially
    weightedDist = exp(-distances/(20^2));
    % normalize
    normDist = weightedDist ./ repmat(sum(weightedDist),length(weightedDist),1);
    % apply interpolation
    interpMap = (templates1005.weights' * normDist)';
    % exclude far off electrode(s)
    matchedTemplate = interpMap(keepElec,:);
end

figure; scatter3(templates1005.chanpos(bestElec,1),templates1005.chanpos(bestElec,2),templates1005.chanpos(bestElec,3),'filled')
hold on; scatter3(elecTrans(keepElec,1),elecTrans(keepElec,2),elecTrans(keepElec,3))
title('Best matching electrodes')


if ~isempty(elecExcludedIndex)
    warning(['Some electrodes have NOT be included in the templates and will need to be removed from your data '...
    'before running templateFit (the list is in electrodesExcludedIndex and electrodesExcludedLabels)'])
end


% re-reference the templates
if refIndex == 0 
    % average reference: substract average across channels
    rerefTemplates = bsxfun(@minus,matchedTemplate, mean(matchedTemplate)); 
else
    % substract the reference channel
    rerefTemplates = bsxfun(@minus,matchedTemplate, matchedTemplate(refIndex,:)); 
end


% output structure
customTemplates.weights = rerefTemplates;
customTemplates.label = chanLabelsData(elecIncludedIndex)';
customTemplates.listROIs = templates1005.listROIs;
customTemplates.reference = refIndex;
customTemplates.matchedLabel = templates1005.label(bestElec(keepElec));
customTemplates.matchedDistance = euclideanDist; % in mm
customTemplates.electrodesIncludedIndex = elecIncludedIndex'; 
customTemplates.electrodesExcludedIndex = elecExcludedIndex';
customTemplates.electrodesExcludedLabel = chanLabelsData(elecExcludedIndex)';


% figure()
% hist(customTemplates.matchedDistance)
% xlabel('Distance in mm')
% ylabel('Number of Electrods')
% title('Distance to HighRes template electrode')


% plot optional topographies 
if eeglabUser && plotMap == 1
    loc = [1:9;10:18]; loc = loc(:);
    mm = round(max(max(abs(customTemplates.weights))),-1);
    figure('position', [200, 200, 1000, 200])
    for roi=1:length(templates1005.listROIs)
        subplot(2,9,loc(roi))
        topoplot(customTemplates.weights(:,roi),channelInfo(elecIncludedIndex),'colormap','jet');
        caxis([-mm mm]);
        title(customTemplates.listROIs(roi));
        axis equal
    end
end  

    