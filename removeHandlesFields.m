function handles = removeHandlesFields (handles)
%function handles = removeHandlesFields (handles)
%
%--------------------------------------------------------------------------
% removeHandlesFields   A simple function to remove some fields from the 
%     handles to be later recalculated 
%
%       INPUT
%         handles:      handles struct with fields
%
%       OUTPUT
%         handles:      clean handles struct
%
%--------------------------------------------------------------------------
%
%     Copyright (C) 2012  Constantino Carlos Reyes-Aldasoro
%
%     This file is part of the PhagoSight package.
%
%     The PhagoSight package is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, version 3 of the License.
%
%     The PhagoSight package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with the PhagoSight package.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
%
% This m-file is part of the PhagoSight package used to analyse fluorescent phagocytes
% as observed through confocal or multiphoton microscopes.  For a comprehensive 
% user manual, please visit:
%
%           http://www.phagosight.org.uk
%
% Please feel welcome to use, adapt or modify the files. If you can improve
% the performance of any other algorithm please contact us so that we can
% update the package accordingly.
%
%--------------------------------------------------------------------------
%
% The authors shall not be liable for any errors or responsibility for the 
% accuracy, completeness, or usefulness of any information, or method in the content, or for any 
% actions taken in reliance thereon.
%
%--------------------------------------------------------------------------
%

if isfield(handles.distanceNetwork, 'forwardRatio')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'forwardRatio'}); 
end
if isfield(handles.distanceNetwork, 'forwardRatioE')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'forwardRatioE'}); 
end
if isfield(handles.distanceNetwork, 'forwardRatioT')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'forwardRatioT'}); 
end
if isfield(handles.distanceNetwork, 'forwardRatioTot')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'forwardRatioTot'}); 
end
if isfield(handles.distanceNetwork, 'inWoundRatio')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'inWoundRatio'}); 
end
if isfield(handles.distanceNetwork, 'inExplorationRatio')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'inExplorationRatio'}); 
end
if isfield(handles.distanceNetwork, 'IdleWoundRatio')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'IdleWoundRatio'}); 
end
if isfield(handles.distanceNetwork,  'LeaveWoundRatio')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'LeaveWoundRatio'}); 
end
if isfield(handles.distanceNetwork, 'LeaveWoundRatio2')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'LeaveWoundRatio2'}); 
end
if isfield(handles.distanceNetwork, 'LeaveWoundRatio3')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'LeaveWoundRatio3'}); 
end

if isfield(handles.distanceNetwork, 'backwardRatioTot')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'backwardRatioTot'}); 
end
if isfield(handles.distanceNetwork, 'numHops')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'numHops'}); 
end
if isfield(handles.distanceNetwork, 'numHopsExploration')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'numHopsExploration'}); 
end
if isfield(handles.distanceNetwork, 'numHopsTranslation')
    handles.distanceNetwork= rmfield(handles.distanceNetwork,{'numHopsTranslation'}); 
end

