% eegplugin_visualareatemplate() - EEGLAB plugin for fitting data using
%                                  visual area templates
%
% Usage:
%   >> eegplugin_visualareatemplate(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Justin Ales, University of St Andrews, 2019


% Copyright (C) 2022 Justin Ales
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

function vers = eegplugin_visualareatemplate(fig, trystrs, catchstrs)

    vers = 'visualareatemplate0.1';
    if nargin < 3
        error('eegplugin_visualareatemplate requires 3 arguments');
    end

    % add folder to path
    % -----------------------
    if ~exist('pop_visualareatemplate')
        p = which('eegplugin_visualareatemplate');
        p = p(1:findstr(p,'eegplugin_visualareatemplate.m')-1);
        addpath([p vers]);
    end
 
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'study');

    % menu callbacks
    % --------------
    comvisualareatemplate = [trystrs.no_check '[STUDY ALLEEG LASTCOM] = pop_visualareatemplate(STUDY, ALLEEG);' catchstrs.new_and_hist];

    % create menus if necessary
    % -------------------------
%     submenu = uimenu( menu, 'Label', 'Visual Area Templates', 'separator', 'on');

    uimenu( menu, 'Label', 'Fit Visual Area Templates', ...
        'CallBack', comvisualareatemplate,...
        'Separator', 'on','userdata', 'startup:on;study:on');

  