function [betaMinNormBest, lambda,residualNorm,solutionNorm,regularizedInverse] = fitTemplates(template , data, plotLcurve, lambda)
% Fit EEG-templates (topographies of functional brain areas) to EEG-data
% INPUTS:
% template = templates for the current EEG montage, either a 2D matrix 
% (electrodes x ROI)or a structure as returned from createCustomTemplates 
% data = average EEG data to fit. For single condition, 2D-matrix with 
% electrodes x time. For multiple conditions, can use a 3-D matrix 
% electrodes x time x condition if all conditions have the same sizes,
% or use a cell array (one cell per condition)
% Time can be replaced by sine & cosine amplitudes when analysing frequencies
% optional plotLcurve: default 0, 1 for plotting L-curve
% optional lambda: regularisation parameter. Will use this value instead of 
% the one computed by the function (use this if the regularisation fails)
% ATTENTION: Make sure that template & data have the same montage and the same reference   
% use createCustomTemplates to fit your EEG montage if needed
% OUTPUTS:
% betaMinNormBest = activity in each ROI for the best regularisation value
% other optional output: 
% lambda = regularisation term, 
% residuals & solutions from the minimum-norm estimates 
% USAGE: roiActivity = fitTemplates(template, data)

addpath('subfunctions')

if (nargin==2), plotLcurve =0; lambda = []; end % do not plot Lcurve & compute best lambda
if (nargin==3), lambda = []; end % plot Lcurve & compute best lambda


% get weights from templates if struct
if isstruct(template)
    template = template.weights;
end


% check if data needs to be reshaped to a 2D matrix
if iscell(data)
    sizeData = cellfun(@size,data,'uni',false);
    data = cell2mat(data);
elseif length(size(data)) == 3 
    nbCond = size(data,3);
    data=reshape(data,[size(data,1),size(data,2)*nbCond]);
end

% check that template and data have the same number of electrodes
% if not, return an error
if size(template,1) ~= size(data,1)
    error(['Mismatch between the number of electrodes in the data and the template. Please check your template.'...
        'Alternatively the data dimensions might be flipped (DIM1 = electrodes, DIM2 = timepoints)'])
end

[u,s,v] = csvd(template);

% use l-curve to find the best regularisation
if isempty(lambda)
    [lambda,~,~,lambdaGridMinNorm] = l_curve_modified(u,s,data,'Tikh',plotLcurve);
    
    if lambda > lambdaGridMinNorm(4) || lambda < lambdaGridMinNorm(end-3)
        warning(['Risk of overfitting. You might want to test another (higher) regularisation value.' ...
            'Try plotting L-curve to determine the corner of the curve using fitTemplates( template , data, 1)'...
            'Then use that corner value in fitTemplates( template , data, 0, ''cornerValue'')'])
    end
end
    
% compute results for the best regularisation
residualNorm = zeros(1,size(data,2));
solutionNorm = zeros(1,size(data,2));
betaMinNormBest = zeros([size(template,2) size(data, 2)]);
for ll = 1:size(data, 2)
    [betaMinNormBest(:, ll),residualNorm(:,ll),solutionNorm(:,ll)] = tikhonov(u, s, v, data(:, ll), lambda);
end


% if multiple conditions, reshape betas to correspond to its input format
if exist('sizeData','var') % input EEG data is cell
    range = cell2mat(arrayfun(@(x) sizeData{x}(2),1:length(sizeData),'uni',false));
    cumul = [0 cumsum(range)];
    betaMinNormBest = arrayfun(@(x) betaMinNormBest(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
    residualNorm = arrayfun(@(x) residualNorm(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
    solutionNorm = arrayfun(@(x) solutionNorm(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
elseif exist('nbCond','var') % input EEG data is 3D mat
    betaMinNormBest = reshape(betaMinNormBest,[size(betaMinNormBest,1),size(betaMinNormBest,2)/nbCond,nbCond]);
    residualNorm = reshape(residualNorm,[size(residualNorm,2)/nbCond,nbCond]);
    solutionNorm = reshape(solutionNorm,[size(solutionNorm,2)/nbCond,nbCond]);
end

%If requested compute Tikhonov regularized inverse matrix
if nargout >=5
    reg_s = diag( s ./ (s.^2 + lambda^2 ));
    regularizedInverse = v * reg_s * u';
end
