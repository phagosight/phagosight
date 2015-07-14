function [handles,numChanges,recordLinks] = linkEndStartTracks (handles,distToLink,sizeGap)
%function [handles,numChanges,recordLinks] = linkEndStartTracks (handles,distToLink,sizeGap)
%
%--------------------------------------------------------------------------
% linkEndStartTracks  analyse tracks as they appear or disappear through time, 
%     try to link an ending track with a starting track, the distance can 
%     be controlled
%
%       INPUT
%         handles:          handles struct including finalNetwork
%         distToLink:       distance (x,y) between node
%         sizeGap:          distance (frames) between nodes
%
%       OUTPUT
%         handles:          updated handles struct
%         numChanges:       number of node linked
%         recordLinks:
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
trackInitNode                       = handles.finalNetwork(1,:);

if ~exist('sizeGap','var')
    sizeGap                         =1;
end
if ~exist('distToLink','var')
    distToLink                      =5;
end

[longestTrack,numTracks]            = size(handles.finalNetwork);
indexFinNodes                       = (sum(handles.finalNetwork>0))+(longestTrack*(0:numTracks-1));
trackFinNode                        = handles.finalNetwork(indexFinNodes);

timesOfFinish                       = sort(unique(handles.nodeNetwork(trackFinNode,5)));

numFinishFr                         = numel(timesOfFinish);
%%
recordLinks                         =[];
numChanges                          =0;
for k=1:numFinishFr

    %%
    currFrame                       =timesOfFinish(k);
    clear tracksFinishHere tracksStartNext positionStart positionFinish
    
    tHasParent                      =  (handles.nodeNetwork(:,8)==0) ;
    tHasChild                       =  (handles.nodeNetwork(:,7)==0);
    tCurrentFrame                   =  (handles.nodeNetwork(:,5)== currFrame   );
    tNextFrame                      =  (handles.nodeNetwork(:,5)==(currFrame+sizeGap));
    tInTrack                        =  (handles.nodeNetwork(:,13)>0);

    tracksFinishHere                = handles.nodeNetwork(tHasParent&tCurrentFrame&tInTrack,13);
    tracksStartNext                 = handles.nodeNetwork(tHasChild&tNextFrame&tInTrack,13);

    NodeTtracksFinishHere           = handles.nodeNetwork(tHasParent&tCurrentFrame&tInTrack,6);
    NodeTracksStartNext             = handles.nodeNetwork(tHasChild&tNextFrame&tInTrack,6);

    positionFinish              = handles.nodeNetwork( trackFinNode( tracksFinishHere),(1:5));
    positionStart               = handles.nodeNetwork( trackInitNode( tracksStartNext),(1:5));

    if ~isempty(positionStart)

        numFinishingT = size(positionFinish,1);
        numStartingT  = size(positionStart,1);
        distanceCube = zeros(numFinishingT,numStartingT);
        distanceCubeZ = zeros(numFinishingT,numStartingT);
        for k2=1:numFinishingT
            distanceCube(k2,:)  = sqrt(sum(((repmat(positionFinish(k2,1:2),[numStartingT 1]) - positionStart(:,1:2)).^2),2));
            distanceCubeZ(k2,:) = sqrt(sum(((repmat(positionFinish(k2,3),[numStartingT 1]) - positionStart(:,3)).^2),2));
        end
        
        % to link any tracks the centroids of the cells must be close in 
        % the XY plane as well as the Z plane!!!
        distanceCube(distanceCubeZ>2)=inf;
        while any(distanceCube(:)<distToLink)
            % Only link if the distances are short.
            % Find which tracks to link and then delete from matrix
            [q1,q2]=find(distanceCube==min(distanceCube(:)));
            q1=q1(1);q2=q2(1);
            numChanges=numChanges+1;
            distanceCube(q1,:)=inf;
            distanceCube(:,q2)=inf;
            
            recordLinks = [recordLinks; [currFrame  tracksFinishHere(q1) tracksStartNext(q2)  trackFinNode(tracksFinishHere(q1))  trackInitNode(tracksStartNext(q2))  ]]; %#ok<AGROW>
            
            % Changes in handles.nodeNetwork, change the track number of the 
            % second track and the finalLabel.
            % Changes in handles.finalNetwork, move nodes from one track to 
            % the end of the other.
            handles.nodeNetwork(handles.nodeNetwork(:,13)==tracksStartNext(q2),13) = tracksFinishHere(q1);

            handles.nodeNetwork(handles.nodeNetwork(:,13)==tracksStartNext(q2),14) = handles.finalLabel( tracksFinishHere(q1));

            %assign parent-child relationship
            handles.nodeNetwork((NodeTtracksFinishHere(q1)) ,8)  = trackInitNode(tracksStartNext(q2));
            handles.nodeNetwork((NodeTracksStartNext(q2)) ,7)  = trackFinNode(tracksFinishHere(q1));

            trackToMove = handles.finalNetwork(:,tracksStartNext(q2));
            trackToMove(trackToMove==0)=[];
            sizeTrackToMove = numel(trackToMove);

            trackToExtend = handles.finalNetwork(:,tracksFinishHere(q1));
            trackToExtend (trackToExtend==0)=[];
            sizeTrackToExtend = numel(trackToExtend);

            %shift track in finalNetwork
            handles.finalNetwork(sizeTrackToExtend+1:sizeTrackToExtend+sizeTrackToMove,tracksFinishHere(q1)) =trackToMove;
            handles.finalNetwork(:,tracksStartNext(q2)) = 0;

            %remove from options to link in the following round
            indexFinNodes([tracksFinishHere(q1) tracksStartNext(q2)])=-1;

        end
    end

end
%%
tracksThatDisappear = find(handles.finalNetwork(1,:)==0);
tracksThatRemain    = find(handles.finalNetwork(1,:)>0);
%%
handles.finalNetwork(:,tracksThatDisappear) =[];
handles.finalLabel(:,tracksThatDisappear)   =[]; 

for counterTR=1:numel(tracksThatRemain)
    handles.nodeNetwork(handles.nodeNetwork(:,13)==tracksThatRemain(counterTR),13) = counterTR;
end

distanceNetwork                 = getDistanceNetForTracks(handles.finalNetwork,handles.nodeNetwork);

handles.distanceNetwork                 = distanceNetwork;
end





