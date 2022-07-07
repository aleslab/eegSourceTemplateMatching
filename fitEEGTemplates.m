function [areaActivity,areaActivityUnscaled,lambda,residualNorm,solutionNorm,regularizedInverse] = fitEEGTemplates(data, templates, plotLcurve, lambda)
%fitEEGTemplates fits EEG-templates (topographies of functional brain areas) to EEG-data
%   A = fitEEGTemplates(data, templates) returns a matrix A containing the
%   scaled (from -1 to 1) contribution of functional brain regions to the EEG signal. 
%   Its format corresponds to the input data (see below).
%
%   DATA is the average (across participants) EEG data. For a single
%   condition, it should be a N-by-T matrix (electrodes x time). For 
%   multiple conditions, it can be a N-by-T-by-C matrix (electrodes x time 
%   x condition) if all conditions have the same size, or a cell array (one
%   cell per condition, each cell in a N-by-T format).
%   Time can be replaced by sine & cosine amplitudes when analysing 
%   frequencies (SSVEP) instead of waveforms.
%
%   TEMPLATES can be a N-by-M matrix (electrodes x regions) or a structure. 
%   The number of electrodes and the reference of the montage should match
%   between the DATA and the TEMPLATES. See the README file and help for
%   createCustomTemplates.
%
%   OPTIONAL INPUT (has to be entered in the following order, use [] to skip one)
%   Useful when the regularisation fails. For more info see README file
%   areaActivity = fitEEGTemplates(myEEGdata, mytemplates [, plotLcurve] [, lambda])
%       'plotLcurve'    - default 0, 1 for plotting L-curve.
%       'lambda'        - regularisation parameter. Will use this value 
%                         instead of the one computed by the function. 
%
%   OPTIONAL OUTPUT
%   [areaActivity, areaActivityUnscaled, lambda,residualNorm,solutionNorm,regularizedInverse] = fitEEGTemplates(data, templates)
%       'areaActivityUnscaled'  - unscaled contribution of each brain regions    
%       'lambda'                - regularisation parameter used in the function. 
%       'residualNorm'          - residuals of the minimum-norm estimates 
%       'solutionNorm'          - solutions of the minimum-norm estimates 
%       'regularizedInverse'    - inverse matrix
%
%   USAGE: 
%       % Returns the contribution of each functional brain regions to myEEGdata
%       areaActivity = fitEEGTemplates(myEEGdata, mytemplates)
%       % Returns the contribution of each functional brain area
%       to myEEGdata using a pre-computed regularisation parameter
%       areaActivity = fitEEGTemplates(myEEGdata, mytemplates, [], 6000)

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

% 
% addpath('subfunctions')

if (nargin==2), plotLcurve =0; lambda = []; end % do not plot Lcurve & compute best lambda
if (nargin==3), lambda = []; end % plot Lcurve & compute best lambda


% get weights from templates if struct
if isstruct(templates)
    templateStruct = templates;
    templates = templates.weights;
end


% check if data needs to be reshaped to a 2D matrix
if iscell(data)
    sizeData = cellfun(@size,data,'uni',false);
    
    %Note tricky use of (:)'.  This is to ensure cell arrays of 
    % [n 1] or [1 n] both concatenate the 2nd dimension of the cells
    % togehter. 
    data = cell2mat(data(:)');
elseif length(size(data)) == 3 
    nbCond = size(data,3);
    data=reshape(data,[size(data,1),size(data,2)*nbCond]);
end

% only keep included electrodes from the template if relevant
if exist('templateStruct','var')
    if ~isempty(templateStruct.electrodesExcludedIndex)
        data = data(templateStruct.electrodesIncludedIndex);
    end
end

% check that templates and data have the same number of electrodes
% if not, return an error
if size(templates,1) ~= size(data,1)
    error(['Mismatch between the number of electrodes in the data and the templates. Please check your templates.'...
        'Alternatively the data dimensions might be flipped (DIM1 = electrodes, DIM2 = timepoints)'])
end

[u,s,v] = csvd(templates);

% use l-curve to find the best regularisation
if isempty(lambda)
    [lambda,~,~,lambdaGridMinNorm] = l_curve_modified(u,s,data,'Tikh',plotLcurve);
    
    if lambda > lambdaGridMinNorm(4) || lambda < lambdaGridMinNorm(end-3)
        warning(['Risk of overfitting. You might want to test another (higher) regularisation value.' ...
            'Try plotting L-curve to determine the corner of the curve using fitEEGTemplates( templates , data, 1)'...
            'Then use that corner value in fitEEGTemplates( templates , data, 0, ''cornerValue'')'])
    end
end
    
% compute results for the best regularisation
residualNorm = zeros(1,size(data,2));
solutionNorm = zeros(1,size(data,2));
areaActivityUnscaled = zeros([size(templates,2) size(data, 2)]);
for ll = 1:size(data, 2)
    [areaActivityUnscaled(:, ll),residualNorm(:,ll),solutionNorm(:,ll)] = tikhonov(u, s, v, data(:, ll), lambda);
end

% normalize betas (contribution) across all ROIs & conditions
areaActivity = areaActivityUnscaled / max(abs(areaActivityUnscaled(:)));

% if multiple conditions, reshape betas to correspond to its input format
if exist('sizeData','var') % input EEG data is cell
    range = cell2mat(arrayfun(@(x) sizeData{x}(2),1:length(sizeData),'uni',false));
    cumul = [0 cumsum(range)];
    areaActivityUnscaled = arrayfun(@(x) areaActivityUnscaled(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
    areaActivity = arrayfun(@(x) areaActivity(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
    residualNorm = arrayfun(@(x) residualNorm(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
    solutionNorm = arrayfun(@(x) solutionNorm(:,cumul(x)+1:cumul(x+1)),1:length(sizeData),'uni',false);
elseif exist('nbCond','var') % input EEG data is 3D mat
    areaActivityUnscaled = reshape(areaActivityUnscaled,[size(areaActivityUnscaled,1),size(areaActivityUnscaled,2)/nbCond,nbCond]);
    areaActivity = reshape(areaActivity,[size(areaActivity,1),size(areaActivity,2)/nbCond,nbCond]);
    residualNorm = reshape(residualNorm,[size(residualNorm,2)/nbCond,nbCond]);
    solutionNorm = reshape(solutionNorm,[size(solutionNorm,2)/nbCond,nbCond]);
end

%If requested compute Tikhonov regularized inverse matrix
if nargout >=5
    reg_s = diag( s ./ (s.^2 + lambda^2 ));
    regularizedInverse = v * reg_s * u';
end
