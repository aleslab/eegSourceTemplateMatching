% pop_visualareatemplate() - Filter data using learned optimal subspaces 
%
% Usage:
%   >> [STUDY, ALLEEG, com] = pop_visualareatemplate(varargin)
%
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%
% Outputs:
%   STUDY        - the input STUDY set with added pre-clustering data for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding 
%                  pre-clustering data (pointers to .mat files that hold cluster measure information).
% Note:

%
% Author: Justin Ales, University of St Andrews, 2022
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

function [STUDY, ALLEEG, com] = pop_visualareatemplate(varargin)

com = '';

if nargin < 1
    help pop_visualareatemplate;
    return
end


if ~ischar(varargin{1}) %initial settings
    if length(varargin) < 2
        error('pop_precomp(): needs both ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    ALLEEG = varargin{2};
    comps = false;
    if nargin > 2
        if strcmpi(varargin{3}, 'components')
            comps = true;
        end
    end
    
    if isempty(ALLEEG)
        error('STUDY contains no datasets');
    end
% if isempty(EEG.data)
%     error('Cannot fit empty dataset.');
% end



% GUI
if nargin < 3
    
  
    
    geometry = {[2 1]};
    geomvert = [1];

    uilist = {{ 'style' 'text' 'string' 'EEG Montage:' } ...
               {'style', 'edit', 'string', 'test','tag','noiseWin'} ...
              };

    [result userdata err structOut] = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Fit Visual Area Templates to Data -- pop_visualareatemplate()', 'helpcom', 'pophelp(''pop_visualareatemplate'')');

    if isempty(result), return; end   

    
else
    
    

    if  nargin <3 || isempty(montageFile)
        error('No montage specified ');
    end
end

% Constants

% Check arguments

%Translate

%Put some checks in here-

betaMinNormBest = fitEEGTemplates(EEG.data(:,:), template);

% EEG.data(:,:) = filterWeights*EEG.data(:,:);
% 
% % In future try to return info about the filter using the ica fields.  
% EEG.icaweights = template;
% EEG.icasphere  = eye(size(template));
% EEG.icawinv    = regularizedInverse;
% EEG.icachansind =  1:EEG.nbchan;

% EEG.icaact = EEG.icaact(windex,:,:);        
        




% fprintf('pop_optimalsubspacefilter() - performing %d point %s filtering.\n', filtorder + 1, filterTypeArray{revfilt + 1, length(edgeArray)})
% fprintf('pop_optimalsubspacefilter() - transition band width: %.4g Hz\n', df)
% fprintf('pop_optimalsubspacefilter() - passband edge(s): %s Hz\n', mat2str(edgeArray))


[STUDY, datavals, times, setinds, cinds] = std_readdata(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(:).labels });

% History string
com = sprintf('%s = pop_optimalsubspacefilter(%s, %s);', inputname(1), inputname(1), vararg2str({method,noiseWin,sigWin,dimensionToFilter}));

end
