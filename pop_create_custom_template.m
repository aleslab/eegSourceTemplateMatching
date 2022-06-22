% pop_create_custom_template() - Create a custom electrode montage template
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

function [STUDY, ALLEEG, com] = pop_create_custom_template(varargin)

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
    
    
    if isempty(ALLEEG)
        error('STUDY contains no datasets');
    end
end


% GUI
if nargin < 3
    
  select_com = ['[inputname, inputpath] = uigetfile2(''*.mat;*.MAT'', ''Choose template file'');'...
              'if inputname ~= 0,' ...
              '   guiind = findobj(''parent'', gcbf, ''tag'', ''templateFile'' );' ...
              '   set( guiind,''string'', [inputpath inputname]);' ...                  
              'end; clear inputname inputpath;'];

    geometry = {[10 30 5]};
    geomvert = [1];

    defaultTemplateFilename = fullfile(STUDY.filepath,'customTemplate.mat')

    uilist = {{ 'style' 'text' 'string' 'Template filename to save:' } ...
               {'style', 'edit', 'string', defaultTemplateFilename,'tag','templateFile'} ...  
               {'style' 'pushbutton' 'string' '...'   'tag' 'browse' 'Callback' select_com} ...
              };



    [result userdata err structOut] = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Fit Visual Area Templates to Data -- pop_visualareatemplate()', 'helpcom', 'pophelp(''pop_visualareatemplate'')');

    if isempty(result), return; end   

    templateFilename = structOut.templateFile
    
else
    
    

    if  nargin <3 || isempty(varargin{3})
        error('No template file specified ');
    end


    templateFilename = varargin{3}
end

% Constants

% Check arguments

%Translate

customTemplate = createCustomTemplates(ALLEEG(1).chanlocs,ALLEEG(1).ref);

save(templateFilename,'-struct','customTemplate');


STUDY.etc.templateFilename = templateFilename;

% History string
com = sprintf('[%s %s] = pop_create_custom_template(%s,%s, %s);', inputname(1),inputname(2),inputname(1),inputname(2), vararg2str({structOut.templateFile}));

end
