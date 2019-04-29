function handles = breakOneTrack(handles,trackToBreak,timeToBreak,woundRegion)
%function handles = breakOneTrack(handles,trackToBreak,timeToBreak,woundRegion)
%
%--------------------------------------------------------------------------
% breakOneTrack  breaks one track into two separate tracks, the first will 
%     be a shorter track and the second will be a new track to be appended  
%     at finalNetwork.
%           Tasks:  1 remove second section of the track,
%                   2 create new track at finalNetwork, finalLabel
%                   2 reassign track parentChild and label in nodeNetwork
%                   1'  recalculate  effectiveTracks2 effectiveDistance x  explorationNetwork, translationNetwork, inWoundNetwork, inWound,
%                   3 remove all calculations from distanceNetwork (recalculate?)
%
%       INPUT
%           handles:        handles structure.
%           trackToBreak:   id number of the track to be broken into two separate tracks.
%           timeToBreak:    time position (frame) at which to break.
%           woundRegion:    if not present then it will be assumed that there is no wound and everything
%                           will be considered TRANSLATION without EXPLORATION
%
%       OUTPUT
%           handles:        new handles with information updated after the break.
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


if (nargin<3)||any(trackToBreak<1)||any(trackToBreak>size(handles.finalNetwork,2))
    help breakOneTrack;
    return;
end

if ~exist('woundRegion','var')
    % if not present then it will be assumed that there is no wound and everything
    % will be considered TRANSLATION without EXPLORATION
    woundRegion                                     = zeros(handles.rows,handles.cols);
end

%%
numTracks = size(handles.finalNetwork,2);
% find the nodes of the track to break
trackIndex = find(handles.finalNetwork (:,trackToBreak));
nodeIndex = handles.finalNetwork(trackIndex,trackToBreak);
trackTimes = handles.nodeNetwork(nodeIndex,5);
positionToBreak = find(trackTimes==timeToBreak);

%[size(handles.distanceNetwork.perHop(...
% positionToBreak:handles.distanceNetwork.numHops(trackToBreak)-1,trackToBreak))]     
% Unsure why in some cases the numHops do not match with the
% size of the finalNetwork, fix it with this
handles.distanceNetwork.numHops = sum(handles.finalNetwork>0);

if isempty(positionToBreak)
    disp('The break time is not within the time of the track');
else
    
    newTrackIndex = nodeIndex(positionToBreak:end);
    lengthNewTrack = numel(newTrackIndex);
    %move node assignment in finalNetwork
    handles.finalNetwork(1:lengthNewTrack,numTracks+1)  = newTrackIndex;
    handles.finalNetwork(positionToBreak:end,trackToBreak)  = 0;
    
    %change parent-child assignment and labels in nodeNetwork
    handles.nodeNetwork(nodeIndex(positionToBreak-1),8)     = 0;
    handles.nodeNetwork(nodeIndex(positionToBreak),7)       = 0;
    handles.nodeNetwork(nodeIndex(positionToBreak),9)       = 0;
    %change label of the new track
    handles.nodeNetwork(nodeIndex(positionToBreak:end),13:14)  = numTracks+1;
    
    %% recalculate distanceNetwork basic parameters:
    handles.distanceNetwork.perHop(1:lengthNewTrack-1,numTracks+1) = ...
        handles.distanceNetwork.perHop(...
            positionToBreak:handles.distanceNetwork.numHops(trackToBreak)-1,trackToBreak);
    handles.distanceNetwork.perHop(positionToBreak:end,trackToBreak) = 0;
    %%
    handles.distanceNetwork.numHops = sum(handles.finalNetwork>0);
    handles.distanceNetwork.totPerTrack = sum(handles.distanceNetwork.perHop);
    handles.distanceNetwork.maxPerTrack = max(handles.distanceNetwork.perHop);
    for k=1:numTracks+1
        handles.distanceNetwork.avPerTrack(k) = ...
            mean(handles.distanceNetwork.perHop(1:handles.distanceNetwork.numHops-1,k));
    end
    handles.distanceNetwork.biasedMean = mean(handles.distanceNetwork.avPerTrack);
    handles.distanceNetwork.biasedStd = std(handles.distanceNetwork.avPerTrack);
    %% Finally calculate the meandering ratio of the tracks if not calculated previously
    distExtBranch(1,numTracks)=0;
    for k=1:numTracks+1
        indexN = handles.finalNetwork(:,k);indexN(indexN==0)=[];
        %calculate now the distance between the extremes of the branch
        distExtBranch(1,k) = sqrt(...
            sum((diff(handles.nodeNetwork(indexN([1 end]),1:3))).^2,2));
    end
    
    % Tortuosity and meanderRatio are inverse of each other
    handles.distanceNetwork.meanderRatio = ...
        distExtBranch./handles.distanceNetwork.totPerTrack;
    handles.distanceNetwork.tortuosity = ...
        handles.distanceNetwork.totPerTrack./distExtBranch;
    
    %%
    if (isfield(handles,'distMaps'))
        handles = rmfield(handles,'distMaps');
    end
    handles = effectiveDistance(handles,woundRegion);
    handles = effectiveTracks(handles,woundRegion);
end
