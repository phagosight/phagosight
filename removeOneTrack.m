function handles = removeOneTrack(handles,trackToRemove,woundRegion)
%function handles = removeOneTrack(handles,trackToRemove)
%function handles = removeOneTrack(handles,trackToRemove,woundRegion)
%
%--------------------------------------------------------------------------
% removeOneTrack  removes a track from the handles
%     Tasks:      1 remove track from finalNetwork, finalLabel
%                 1'  recalculate  effectiveTracks2 effectiveDistance x  explorationNetwork, translationNetwork, inWoundNetwork, inWound,
%                 2 remove track assignment in nodeNetwork
%                 3 remove all calculations from distanceNetwork (recalculate?)
%
%       INPUT
%         handles:          handles struct
%         trackToRemove:    label (id number) of the track to remove
%         woundRegion:      wound region matrix (to recalculate metrics)
%
%       OUTPUT
%         handles:          updated handles
           
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
% accuracy, completeness, or usefulness of any information, or method in the
% content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

if (nargin<2)||any(trackToRemove<1)||any(trackToRemove>size(handles.finalNetwork,2))
    % help removeOneTrack;
    fprintf('%s: Error with trackToRemove: %d \n', mfilename, trackToRemove);
    return;
end

if ~exist('woundRegion','var')
    % if not present then it will be assumed that there is no wound and everything
    % will be considered TRANSLATION without EXPLORATION
    woundRegion                                     = zeros(handles.rows,handles.cols);
end

% find the nodes of the track to remove
trackIndex                                  = find(handles.finalNetwork (:,trackToRemove));
nodeIndex                                   = handles.finalNetwork(trackIndex,trackToRemove);

handles.finalNetwork(:,trackToRemove)           = [];


handles.finalLabel                          = neutrophilOverlap(handles.nodeNetwork,handles.finalNetwork,handles);

handles.nodeNetwork(:,[9 13:14])                = 0;

for counterTracks=1:size(handles.finalNetwork,2)
    trackIndex                                  = find(handles.finalNetwork (:,counterTracks));
    nodeIndex                                   = handles.finalNetwork(trackIndex,counterTracks); %#ok<FNDSB>
    handles.nodeNetwork(nodeIndex,13)           = counterTracks; %#ok<FNDSB>
    handles.nodeNetwork(nodeIndex,14)           = handles.finalLabel(counterTracks); %#ok<FNDSB>
    handles.nodeNetwork(nodeIndex(2:end),9)     = sqrt(sum(diff(handles.nodeNetwork(nodeIndex,1:3)).^2,2));
end


try handles.distanceNetwork.perHop(:,trackToRemove)                     = []; catch ;    qqq=0;end
try handles.distanceNetwork.totPerTrack(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.maxPerTrack(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.avPerTrack(:,trackToRemove)                 = []; catch ;    qqq=0;end
try handles.distanceNetwork.numHopsInWound(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.forwardRatioTot2(:,trackToRemove)           = []; catch ;    qqq=0;end
try handles.distanceNetwork.staticRatio(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.returnRatio(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.returnRatio2(:,trackToRemove)               = []; catch ;    qqq=0;end
try handles.distanceNetwork.absVelocity(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.absVelocityExp(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.absVelocityTra(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.oriVelocity(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.oriVelocityExp(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.oriVelocityTra(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.latVelocity(:,trackToRemove)                = []; catch ;    qqq=0;end
try handles.distanceNetwork.latVelocityExp(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.latVelocityTra(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.crossWoundBack(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.crossWoundBackRatio(:,trackToRemove)        = []; catch ;    qqq=0;end
try handles.distanceNetwork.numHops(:,trackToRemove)                    = []; catch ;    qqq=0;end
try handles.distanceNetwork.numHopsExploration(:,trackToRemove)         = []; catch ;    qqq=0;end
try handles.distanceNetwork.numHopsTranslation(:,trackToRemove)         = []; catch ;    qqq=0;end
try handles.distanceNetwork.forwardRatioTot(:,trackToRemove)            = []; catch ;    qqq=0;end
try handles.distanceNetwork.forwardRatioE(:,trackToRemove)              = []; catch ;    qqq=0;end
try handles.distanceNetwork.forwardRatioT(:,trackToRemove)              = []; catch ;    qqq=0;end
try handles.distanceNetwork.inExplorationRatio(:,trackToRemove)         = []; catch ;    qqq=0;end
try handles.distanceNetwork.inWoundRatio(:,trackToRemove)               = []; catch ;    qqq=0;end
try handles.distanceNetwork.IdleWoundRatio(:,trackToRemove)             = []; catch ;    qqq=0;end
try handles.distanceNetwork.backwardRatioTot(:,trackToRemove)           = []; catch ;    qqq=0;end
try handles.distanceNetwork.LeaveWoundRatio(:,trackToRemove)            = []; catch ;    qqq=0;end
try handles.distanceNetwork.LeaveWoundRatio2(:,trackToRemove)           = []; catch ;    qqq=0;end
try handles.distanceNetwork.LeaveWoundRatio3(:,trackToRemove)           = []; catch ;    qqq=0;end

try handles.explorationNetwork(:,trackToRemove)                         = []; catch ;    qqq=0;end
try handles.translationNetwork(:,trackToRemove)                         = []; catch ;    qqq=0;end
try handles.inWoundNetwork(:,trackToRemove)                             = []; catch ;    qqq=0;end
try handles.inWound(:,trackToRemove)                                    = []; catch ;    qqq=0;end

handles.distanceNetwork.biasedMean                                  = mean(handles.distanceNetwork.avPerTrack);
handles.distanceNetwork.biasedStd                                   = std(handles.distanceNetwork.avPerTrack);

if exist('woundRegion','var')
    if isempty(woundRegion)
        %handles                                                         = effectiveTracks(handles);
        %handles                                                     = rmfield(handles,'distMaps');
    else
        if (isfield(handles,'distMaps'))
            handles                                                     = rmfield(handles,'distMaps');
        end
        handles                                                         = effectiveDistance(handles,woundRegion);
        handles                                                         = effectiveTracks(handles,woundRegion);
    end
else
    if (isfield(handles,'distMaps'))
        handles.distMaps.absDistPerHop(:,trackToRemove)                 = [];
        handles.distMaps.anglePerHop(:,trackToRemove)                   = [];
        handles.distMaps.oriDistPerHop(:,trackToRemove)                 = [];
        handles.distMaps.furthestPoint(:,trackToRemove)                 = [];
        handles.distMaps.furthestPoint2(:,trackToRemove)                = [];
        handles.distMaps.latDistPerHop(:,trackToRemove)                 = [];
        handles.distMaps.furthestPoint3(:,trackToRemove)                = [];
        handles.distMaps.effDistPerHop(:,trackToRemove)                 = [];
        handles.distMaps.initialPoint3(:,trackToRemove)                 = [];
        handles.distMaps.initialPoint2(:,trackToRemove)                 = [];
        handles.distMaps.activationPoint3(:,trackToRemove)              = [];
        handles.distMaps.activationTime3(:,trackToRemove)               = [];
        handles.distMaps.activationDistTW(:,trackToRemove)              = [];
        handles.distMaps.enterWoundPoint3(:,trackToRemove)              = [];
        handles.distMaps.enterWoundTime3(:,trackToRemove)               = [];
    end
end


