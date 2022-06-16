function customTemplates = createCustomTemplates(channelInfo,varargin)
% function creating EEG topographic templates according to the EEG montage  
% with which the dataset was recorded
% Plot the topographies of the customTemplates for each ROI to check that
% they are similar to what is in the paper. Done automatically for EEGlab
% data. If it does not seem right, check the coordinate system of
% your data (see further help in the function)
% example EEGlab: customTemplates = createCustomTemplates(EEG.chanlocs)
% example Fieldtrip: customTemplates = createCustomTemplates(cfg.elec)
%
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
%
% output is a structure containing:
%   customTemplates.weights = matrix of EEG response for N electrodes x 18 ROI-templates 
%   customTemplates.listROIs = cell-array of the 18 ROI names (for L-left and R-right separately)
%   customTemplates.chanLabels = cell-array of length N with the label of
%   each electrode (same as the dataset)
%   customTemplates.reference = reference used for the templates (0=average)
%   customTemplates.matchedLabels = best matching labels from the templates
%   customTemplates.matchedDistances = euclidean distances (in mm) between
%   electrodes comparing user & standard EEG montage  
%   customTemplates.electrodesIncludedIndex = electrodes' indexes
%   corresponding to the chanLabels included in the custom templates
%   customTemplates.electrodesExcludedIndex = indexes of excluded
%   electrodes from the user's montage (due to missing location information 
%   or no satisfactory matching)
%   customTemplates.electrodesExcludedLabels = labels of excluded
%   electrodes
% Optional input (has to be in the following order, [] to skip one):
% 1st: Reference of the EEG montage. Default is read from channelInfo for 
% EEGlab. If not found, ask user or can be entered here. The reference 
% should be the channel number or 0 for average reference. 
% usage: customTemplates = createCustomTemplates(EEG.chanlocs,0)
% 2nd: plot the templates 1 (default) or not 0. Plotting only works for 
% EEG.chanlocs input
% 3rd: Coordinate system: 'ALS' or 'RAS'. By default the coordinate system  
% is ALS for EEGlab users, RAS otherwise but can be specified by the user. 
% (see below for more details. Other coordinate system can potentially be 
% added in the program) 
% usage: customTemplates = createCustomTemplates(EEG.chanlocs,[],[],'RAS')
% 4th: interpolation = default 0, set to 1 to use interpolation method
% instead of closest electrode to fit the montage
%
% To do the fitting, need 3-D cartesian coordinates of the electrodes
% corresponding to the EEG montage that is used AND at least 8 channel 
% labels to match. These can be extracted from the dataset or entered 
% manually when prompted if no label is found.
% Also need the reference of the montage
%
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

addpath('subfunctions')

% 4 maximum optional inputs
if length(varargin)>4
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
    disp('%d electrodes do not have location data',sum(~tf))
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
while length(refIndex)~=1 || ~ismember(refIndex,0:length(chanLabelsData)) % include 0 for average reference
    msgError = errordlg('Please enter a valid number for the EEG reference','Input error','modal');
    uiwait(msgError);
    prompt = {'Enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create custom templates';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end
    
% load standard 10-05 templates 
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


figure()
hist(customTemplates.matchedDistance)
xlabel('Distance in mm')
ylabel('Number of Electrods')
title('Distance to HighRes template electrode')


% plot optional topographies 
if eeglabUser && plotMap == 1
    loc = [1:9;10:18]; loc = loc(:);
    mm = round(max(max(abs(customTemplates.data))),-1);
    figure('position', [200, 200, 1000, 200])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(customTemplates.data(:,roi),channelInfo,'colormap','jet');
        caxis([-mm mm]);
        title(listROIs(roi));
        axis equal
    end
end  

    