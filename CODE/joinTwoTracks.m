function handles = joinTwoTracks(handles,trackOne,trackTwo,...
                                 woundRegion,limitDistance,limitTime)
%function handles = joinTwoTracks(handles,trackOne,trackTwo,woundRegion)
%
%--------------------------------------------------------------------------
% joinTwoTracks  two tracks will be merged into a single one in the present 
%     handles.
%       Tasks:
%           1 pass nodes from trackTwo to trackOne
%           2 remove trackTwo from finalNetwork
%           3 recalculate  effectiveTracks effectiveDistance x
%           explorationNetwork,translationNetwork, inWoundNetwork, inWound,
%           4 reassign track assignment in nodeNetwork
%           5 remove all calculations from distanceNetwork (recalculate?)
%
%       INPUT
%         handles:              handles struct
%         trackOne,trackTwo:    id numbers of the tracks to be merged
%         woundRegion:          wound region (to recalculate metrics)
%
%       OUTPUT
%         handles:              new handles struct with the tracks merged
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
% This m-file is part of the PhagoSight package used to analyse
% fluorescent phagocytes as observed through confocal or
% multiphoton microscopes.  For a comprehensive user manual, please visit:
%
%           <http://www.phagosight.org.uk>
%
% Please feel welcome to use, adapt or modify the files. If you can improve
% the performance of any other algorithm please contact us so that we can
% update the package accordingly.
%
%--------------------------------------------------------------------------
%
% The authors shall not be liable for any errors or responsibility for the 
% accuracy, completeness, or usefulness of any information, or
% method in the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

%%
numTracks = size(handles.finalNetwork,2);

if (nargin<3)||any(trackOne<1)||any(trackOne>numTracks)||...
        any(trackTwo<1)||any(trackTwo>numTracks)
    help joinTwoTracks;
    return;
end
if ~exist('woundRegion','var')
    % if not present then it will be assumed that there is no wound and everything
    % will be considered TRANSLATION without EXPLORATION
    woundRegion = zeros(handles.rows,handles.cols);
end
if ~exist('limitTime','var')
    limitTime = 16;
end
%%
%Condition for now, they have to be consecutive (t, t+1) or overlapping (t,t)
if trackTwo<trackOne
    trackTemp = trackOne;
    trackOne  = trackTwo;
    trackTwo  = trackTemp;
end

trackIndex1 = find(handles.finalNetwork (:,trackOne));
nodeIndex1 = handles.finalNetwork(trackIndex1,trackOne);
timeIndex1 = handles.nodeNetwork(nodeIndex1,5);
numNodesTrack1 = numel(trackIndex1);

trackIndex2 = find(handles.finalNetwork (:,trackTwo));
nodeIndex2 = handles.finalNetwork(trackIndex2,trackTwo);
timeIndex2 = handles.nodeNetwork(nodeIndex2,5);
numNodesTrack2 = numel(trackIndex2);

%it may be the case that a track with a higher label is before one with a
%larger label, re-order in that case
if timeIndex1(1)>=timeIndex2(end)
    trackIndexTemp = trackIndex1;
    nodeIndexTemp = nodeIndex1;
    timeIndexTemp = timeIndex1;
    
    
    trackIndex1 = trackIndex2;
    nodeIndex1 = nodeIndex2;
    timeIndex1 = timeIndex2;
    trackIndex2 = trackIndexTemp;
    nodeIndex2 = nodeIndexTemp;
    timeIndex2 = timeIndexTemp;
    
    trackTemp = trackOne;
    trackOne  = trackTwo;
    trackTwo  = trackTemp;
end
joinTheTracks =0;

%%
%Calculating how time overlap, how many frame have the tracks in common
timeDistance = timeIndex2(1) - timeIndex1(end);

% Only joint if the tracks do not overlap completely (First if) AND
% Only joint if there is a overlap in the range [-4 16]  (second if)

%if (numel(union(timeIndex2,timeIndex1))>max(2+numNodesTrack1,2+numNodesTrack2))
if (numel(intersect(timeIndex2,timeIndex1))<min(numNodesTrack1,numNodesTrack2))
    
    if ((timeDistance > -5) && (timeDistance < limitTime))
        
        if (timeDistance < 1)
            idxRef = 2 + (-1 * timeDistance);
            try
                distance = (handles.nodeNetwork(nodeIndex1(end),1) ...
                            - handles.nodeNetwork(nodeIndex2(...
                                min(numel(nodeIndex2),idxRef)),1)).^2;
            catch
                qqq = 1;
            end
            
            distance = distance + (handles.nodeNetwork(nodeIndex1(end),2) - ...
                                   handles.nodeNetwork(nodeIndex2(...
                                       min(numel(nodeIndex2),idxRef)),2)).^2;
            distance = sqrt(distance);
            pctTime = 1;
        else
            distance = (handles.nodeNetwork(nodeIndex1(end),1) - ...
                        handles.nodeNetwork(nodeIndex2(1),1)).^2;
            distance = distance + (handles.nodeNetwork(nodeIndex1(end),2) - ...
                                   handles.nodeNetwork(nodeIndex2(1),2)).^2;
            distance = sqrt(distance);
            pctTime = timeDistance;
        end
        
        pctTime = pctTime/handles.numFrames;
        pctPxl = 9/10 * pctTime + .1;
        if ~exist('limitDistance','var')
            limitDistance = pctPxl * max(handles.rows,handles.cols);
        end
        %fprintf('%d frames -> %f px Limit | %f
        %px\n',timeDistance,...
        % limitDistance,distance);
        
        if (distance < limitDistance)
            joinTheTracks = 1;
            if (timeDistance < 1)
                idxRef = 2 + (-1 * timeDistance);
                %index of first track in final network
                idxBeg = 1 + trackIndex1(end);
                idxEnd = trackIndex1(end) + trackIndex2(end) + (timeDistance - 1);
                %Nodes in trackTwo (nodeIndex2) are moved to the trackOne in
                %the finalNetwork
                handles.finalNetwork(idxBeg:idxEnd,trackOne) = ...
                    nodeIndex2(idxRef:end);
                %assign parent-child
                handles.nodeNetwork(nodeIndex1(end),8) = nodeIndex2(idxRef);
                handles.nodeNetwork(nodeIndex2(idxRef),7) = nodeIndex1(end);
                handles.nodeNetwork(nodeIndex2(1:(idxRef-1)),7:8) = 0;
                
                %change track assignment
                handles.nodeNetwork(nodeIndex2(idxRef:end),13) = ...
                    handles.nodeNetwork(nodeIndex1(end),13);
                handles.nodeNetwork(nodeIndex2(idxRef:end),14) = ...
                    handles.nodeNetwork(nodeIndex1(end),14);
                handles.nodeNetwork(nodeIndex2(1:(idxRef-1)),13:14) = 0;
                
                %recalculate parameters in distanceNetwork only for new track
                trackIndex3 = [trackIndex1;trackIndex1(end)+...
                               trackIndex2(1:end-(idxRef-1))];
                numNodesBranches = trackIndex3(end);
                indexN = handles.finalNetwork(trackIndex3,trackOne);
                
            else
                %Nodes in trackTwo (nodeIndex2) are moved to the trackOne in
                %the finalNetwork
                handles.finalNetwork(1+trackIndex1(end):trackIndex1(end)+...
                                     trackIndex2(end),trackOne) = nodeIndex2;
                
                %assign parent-child
                handles.nodeNetwork(nodeIndex1(end),8) = nodeIndex2(1);
                handles.nodeNetwork(nodeIndex2(1),7) = nodeIndex1(end);
                
                %change track assignment
                handles.nodeNetwork(nodeIndex2,13) = ...
                    handles.nodeNetwork(nodeIndex1(end),13);
                handles.nodeNetwork(nodeIndex2,14) = ...
                    handles.nodeNetwork(nodeIndex1(end),14);
                %recalculate parameters in distanceNetwork only for new track
                trackIndex3 = [trackIndex1;trackIndex1(end)+trackIndex2];
                numNodesBranches = trackIndex3(end);
                indexN = handles.finalNetwork(trackIndex3,trackOne);
            end
            
        else
            STR = ['The spatial distance (' num2str(distance) ...
                   ' pixels) between the tracks is too large to join them. ' ...
                   'For a temporal separation of ' num2str(timeDistance) ...
                   ' frames, the distance limit is ' ...
                   num2str(limitDistance) ' pixels.'];
            disp('ERROR');
            disp(STR);
        end
        
    else
        STR = ['The temporal distance (' num2str(timeDistance) ...
               ' frames) between the tracks is too large to join them. '...
               'It is only possible to joint tracks with a temporal ' ...
               'distance between -4 and 15 frames.'];
        disp('ERROR');
        disp(STR);
    end
else
    STR = ['The tracks overlap completely in time' ...
           'It is only possible to joint tracks that continue each ' ...
           'other.'];
    disp('ERROR');
    disp(STR);
end

%%

if joinTheTracks==1
    
    handles.distanceNetwork.perHop(trackIndex3(1:end-1),trackOne) = ...
        sqrt(sum((diff(handles.nodeNetwork(indexN,1:3))).^2,2)); 
    handles.distanceNetwork.numHops(1,trackOne) = numNodesBranches;
    handles.distanceNetwork.totPerTrack(1,trackOne) = ...
        sum(handles.distanceNetwork.perHop(trackIndex3(1:end-1),trackOne));
    handles.distanceNetwork.maxPerTrack(1,trackOne) = ...
        max(handles.distanceNetwork.perHop(trackIndex3(1:end-1),trackOne));
    handles.distanceNetwork.avPerTrack(1,trackOne) =...
        handles.distanceNetwork.totPerTrack(1,trackOne)./(numNodesBranches-1);
    
    handles = removeOneTrack(handles,trackTwo,woundRegion);
    
    if isfield(handles,'distMaps')
        handles = rmfield(handles,'distMaps');
        if exist('woundRegion','var')
            if isempty(woundRegion)
                %handles = effectiveDistance(handles);
            else
                handles = effectiveDistance(handles,woundRegion);
                handles = effectiveTracks(handles,woundRegion);
            end
        else
            handles = effectiveDistance(handles);
        end
    end
end
   
