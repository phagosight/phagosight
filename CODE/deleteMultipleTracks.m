function handles = deleteMultipleTracks(handles,woundRegion,minDistanceTrack)
%function handles = deleteMultipleTracks(handles,woundRegion,minDistanceTrack)
%
%
%--------------------------------------------------------------------------
% deleteMultipleTracks delete all the tracks that are below a certain distance
%      and calls removeOneTrack for every candidate
%       Tasks:
%           1 pass nodes from trackTwo to trackOne
%           2 remove trackTwo from finalNetwork
%           3 recalculate  effectiveTracks effectiveDistance x  explorationNetwork, translationNetwork, inWoundNetwork, inWound,
%           4 reassign track assignment in nodeNetwork
%           5 remove all calculations from distanceNetwork (recalculate?)
%
%       INPUT
%         handles:                  handles struct
%         woundRegion:              wound region (to recalculate metrics)
%         maxDistanceConsidered:    above the maximum distance the pair is discarded
%         minDistanceConsidered:    below the minimum distance the pair is discarded
%         minSeparationTwin:        in case one track may be assigned as child to two
%                                   candidates, only assign if one is closer to the
%                                   child by this distance (default is to discard
%                                   both options and be done manually)
%
%       OUTPUT
%         handles:                  new handles struct with the tracks merged
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

%%

numTracks = size(handles.finalNetwork,2);

if ~exist('woundRegion','var')
    woundRegion = zeros(handles.rows,handles.cols);
    woundRegion (:,ceil(handles.cols/4):handles.cols) = 1;
end
%
%listTracks = 1:numTracks;

%%
recordTrackDistances =[];
for counterTrack = numTracks:-1:1  %39;
    
    currTrack                   = handles.nodeNetwork(handles.nodeNetwork(:,13)==counterTrack,1:5);
    % Only consider tracks that are shorter than a certain size.
    if size(currTrack,1) <minDistanceTrack
        handles                 = removeOneTrack(handles,counterTrack,woundRegion);

        
    end
    
end
end
