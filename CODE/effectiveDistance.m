function handles=effectiveDistance(handles, woundRegion)
%function handles = effectiveDistance(handles)
%function handles = effectiveDistance(handles, woundRegion)
%
%--------------------------------------------------------------------------
% effectiveDistance  calculates the effective distance that is traversed by
%     a tracked object the orientation of the track is calculated from the 
%     first and last points and all points are rotated to project the 
%     movement to the orientation. Then, the distance is calculated.
%     If the woundRegion is provided, it then calculates whether a neutrophil 
%     has entered the region and all tracks are divided into two states:
%     Translation all movement before it reaches the wound for the first time
%     Exploration all movement after it has reached the wound.
%
%
%     The calculation of effectiveDistance adds handles.translationNetwork 
%     (those movements before the neutrophil reaches the wound), 
%     handles.explorationNetwork (those movements after the neutrophil 
%     reaches the wound) and handles.inWoundNetwork (a register if the 
%     neutrophil is inside or outside the wound). 
%     In handles.distanceNetwork, the following fields are created: 
%
%     perHopT,             perHopE,             perHopOrientedT,
%     perHopEffectiveT,    perHopOrientedE,     perHopEffectiveE,
%     numHopsTranslation,  numHopsExploration,  forwardRatio,
%     forwardRatioE,       forwardRatioTot,     inWoundRatio,
%     inExplorationRatio,  IdleWoundRatio,      LeaveWoundRatio,
%     LeaveWoundRatio2
%          
%   
%       INPUT
%         handles:      handles structure.
%         woundRegion:  matrix indicating the wound region
%
%       OUTPUT
%         handles:      handles structure updated with the fields described
%                       above
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


numTracks = size(handles.distanceNetwork.perHop,2);

% A health check, if there are tracks with one cell (no hops at all) remove.

tracksToRemove                      = (handles.distanceNetwork.numHops==0);

while any(tracksToRemove>0)
   % if(exist('woundRegion','var'))
   %     handles                     = removeOneTrack(handles,find(tracksToRemove,1),woundRegion);
   % else
        handles                     = removeOneTrack(handles,find(tracksToRemove,1));
   % end
   numTracks = size(handles.distanceNetwork.perHop,2);
    tracksToRemove                  = (handles.distanceNetwork.numHops==0);
end

%If wound region is provided, split tracks into two regions Translation and Exploration (of the wound)
if exist('woundRegion','var')

    [H]                                                                 = hough(woundRegion);
    peaks                                                               = houghpeaks(H, 1);
    woundAngle                                                          = pi*(0+(180-peaks(2)))/180;
    %%
    if isfield(handles.distanceNetwork,'perHopEffectiveE')
        handles.distanceNetwork= rmfield(handles.distanceNetwork,{'perHopEffectiveT','perHopEffectiveE','perHopT','perHopE','perHopOrientedT','perHopOrientedE'});
        handles.distanceNetwork= rmfield(handles.distanceNetwork,{'inWoundRatio','inExplorationRatio','IdleWoundRatio','LeaveWoundRatio','LeaveWoundRatio2'});
        handles.distanceNetwork= rmfield(handles.distanceNetwork,{'numHopsExploration','numHopsTranslation'});
        
    end
    %%
    handles.distanceNetwork.perHopT(1,numTracks)                        = 0;
    handles.distanceNetwork.perHopE(1,numTracks)                        = 0;
    handles.distanceNetwork.perHopOrientedT(1,numTracks)                = 0;
    handles.distanceNetwork.perHopEffectiveT(1,numTracks)               = 0;
    handles.explorationNetwork(1,numTracks)                             = 0;
    handles.translationNetwork(1,numTracks)                             = 0;
    handles.distanceNetwork.perHopOrientedE(1,numTracks)                = 0;
    handles.distanceNetwork.perHopEffectiveE(1,numTracks)               = 0;
    %First detect which regions of the track are before first visit to the wound
    for currTrack = 1:numTracks
        %find how long is the track
        lastNodeT                                                                   = handles.distanceNetwork.numHops(currTrack);
       
        %to determine if a point of a track is within the wound, get the index for the wound and the
        %index for the track nodes and compare them with ismember
        indexWound                                                                  = find (woundRegion);
        tempIndex                                                                   = sub2ind([handles.rows handles.cols],round(handles.nodeNetwork(handles.finalNetwork(1:lastNodeT,currTrack),1)),round(handles.nodeNetwork(handles.finalNetwork(1:lastNodeT,currTrack),2)));
        inWound                                                                     = (ismember(tempIndex,indexWound));

        % inWound has a logical index: 0 if the node is not in the wound, 1 otherwise, determine if the
        % track ever starts exploring
        startExploration                                                            = (find(inWound,1,'first')-1);


        %new fields:    translationNetwork/explorationNetwork same structure as finalNetwork
        %               numHopsTranslation/exploration          number of nodes
        %               perHopT/E                               The distance of each hop

        if isempty(startExploration)
            %no exploration at all, only translation, fill in details in new fields
            inTranslation                                                               = 1:lastNodeT;
            handles.translationNetwork(inTranslation,currTrack)                         = handles.finalNetwork(inTranslation,currTrack);
            handles.distanceNetwork.numHopsTranslation(currTrack)                       = lastNodeT-1;
            handles.distanceNetwork.numHopsExploration(currTrack)                   	= 0;
            
            %keep the index inWound for future reference
            handles.inWoundNetwork(inTranslation,currTrack)                             = inWound (inTranslation);
            initFinCoord                                                                = handles.nodeNetwork(handles.finalNetwork((inTranslation),currTrack),(1:2));

            tempAngle                                                                   = angle((initFinCoord(end,1)-initFinCoord(1,1))+i*(initFinCoord(end,2)-initFinCoord(1,2)));
            tempXY                                                                      = ([cos(tempAngle) sin(tempAngle);-sin(tempAngle) cos(tempAngle)]*initFinCoord(:,[1 2])')';
            tempD                                                                       = diff(tempXY); %#ok<UDIM>
            
            handles.distanceNetwork.perHopT(inTranslation(1:end-1),currTrack)           = sqrt(sum(diff(initFinCoord).^2,2));
            handles.distanceNetwork.perHopOrientedT(inTranslation(1:end-1),currTrack)   = tempD(:,1);
            handles.distanceNetwork.perHopEffectiveT(inTranslation(1:end-1),currTrack)  = handles.distanceNetwork.perHopOrientedT(inTranslation(1:end-1),currTrack)./handles.distanceNetwork.perHopT(inTranslation(1:end-1),currTrack);

        elseif startExploration==0
            %no translation at all, only exploration
            inExploration                                                           = 1+startExploration:lastNodeT;
            handles.explorationNetwork(inExploration-startExploration,currTrack)    = handles.finalNetwork(inExploration,currTrack);
            handles.distanceNetwork.numHopsTranslation(currTrack)                   = startExploration;
            handles.distanceNetwork.numHopsExploration(currTrack)                   = lastNodeT - 1;
            
             
            %keep the index inWound for future reference
            handles.inWoundNetwork(inExploration-startExploration,currTrack)        = inWound (inExploration);

            %calculate effective distances, second for the exploration, from limit to the wound to end
            initFinCoord                                                                = handles.nodeNetwork(handles.finalNetwork((inExploration),currTrack),(1:2));

            tempXY                                                                      = ([cos(woundAngle) sin(woundAngle);-sin(woundAngle) cos(woundAngle)]*initFinCoord(:,[1 2])')';
            tempD                                                                       = diff(tempXY); %#ok<UDIM>
            
            handles.distanceNetwork.perHopE(inExploration(1:end-1),currTrack)           = sqrt(sum(diff(initFinCoord).^2,2));
            handles.distanceNetwork.perHopOrientedE(inExploration(1:end-1),currTrack)   = tempD(:,1);
            handles.distanceNetwork.perHopEffectiveE(inExploration(1:end-1),currTrack)  = handles.distanceNetwork.perHopOrientedE(inExploration(1:end-1),currTrack)./handles.distanceNetwork.perHopE(inExploration(1:end-1),currTrack);
   
            
        else
            %Divide the track into translation and exploration regions
            inTranslation                                                           = 1:startExploration;
            inExploration                                                           = 1+startExploration:lastNodeT;
            handles.translationNetwork(inTranslation,currTrack)                     = handles.finalNetwork(inTranslation,currTrack);
            handles.explorationNetwork(inExploration-startExploration,currTrack)    = handles.finalNetwork(inExploration,currTrack);
            handles.distanceNetwork.numHopsTranslation(currTrack)                   = startExploration;
            handles.distanceNetwork.numHopsExploration(currTrack)                   = lastNodeT - startExploration;
            
            %keep the index inWound for future reference
            handles.inWoundNetwork(inExploration-startExploration,currTrack)        = inWound (inExploration);

            
            %calculate effective distances, first for the translation, from 1 to limit to the wound
            initFinCoord                                                                = handles.nodeNetwork(handles.finalNetwork((inTranslation),currTrack),(1:2));
            
            tempAngle                                                                   = angle((initFinCoord(end,1)-initFinCoord(1,1))+i*(initFinCoord(end,2)-initFinCoord(1,2)));
            
            tempXY                                                                      = ([cos(tempAngle) sin(tempAngle);-sin(tempAngle) cos(tempAngle)]*initFinCoord(:,[1 2])')';
            tempD                                                                       = diff(tempXY); %#ok<UDIM>
            
            handles.distanceNetwork.perHopT(inTranslation(1:end-1),currTrack)           = sqrt(sum(diff(initFinCoord).^2,2));
            handles.distanceNetwork.perHopOrientedT(inTranslation(1:end-1),currTrack)   = tempD(:,1);
            handles.distanceNetwork.perHopEffectiveT(inTranslation(1:end-1),currTrack)  = handles.distanceNetwork.perHopOrientedT(inTranslation(1:end-1),currTrack)./handles.distanceNetwork.perHopT(inTranslation(1:end-1),currTrack);
 
            %calculate effective distances, second for the exploration, from limit to the wound to end
            initFinCoord                                                                = handles.nodeNetwork(handles.finalNetwork((inExploration),currTrack),(1:2));

            tempXY                                                                      = ([cos(woundAngle) sin(woundAngle);-sin(woundAngle) cos(woundAngle)]*initFinCoord(:,[1 2])')';
            tempD                                                                       = diff(tempXY); %#ok<UDIM>
            
            handles.distanceNetwork.perHopE(inExploration(1:end-1)-startExploration,currTrack)           = sqrt(sum(diff(initFinCoord).^2,2));
            handles.distanceNetwork.perHopOrientedE(inExploration(1:end-1)-startExploration,currTrack)   = tempD(:,1);
            handles.distanceNetwork.perHopEffectiveE(inExploration(1:end-1)-startExploration,currTrack)  = handles.distanceNetwork.perHopOrientedE(inExploration(1:end-1)-startExploration,currTrack)./handles.distanceNetwork.perHopE(inExploration(1:end-1)-startExploration,currTrack);
              
            
        end

    end


    %% Forward ratio in translation
    %calculate how forward-idle-backwards balance for each track on the TRANSLATION
    tempDenom                               = handles.distanceNetwork.numHopsTranslation;
    tempDenom(tempDenom==0)                 = inf;
    handles.distanceNetwork.forwardRatio    = sum(handles.distanceNetwork.perHopEffectiveT>0.3)./tempDenom;
    %% Forward ratio in exploration
    %calculate how forward-idle-backwards balance for each track on the Exploration
    tempDenom                               = handles.distanceNetwork.numHopsExploration;
    tempDenom(tempDenom==0)                 = inf;
    handles.distanceNetwork.forwardRatioE    = sum(handles.distanceNetwork.perHopEffectiveE>0.3)./tempDenom;
    %% Forward ratio in total
    %calculate how forward-idle-backwards balance for each track on the Exploration
    tempDenom                               = handles.distanceNetwork.numHops;
    tempDenom(tempDenom==0)                 = inf;
    handles.distanceNetwork.forwardRatioTot    = (sum(handles.distanceNetwork.perHopEffectiveE>0.3)+sum(handles.distanceNetwork.perHopEffectiveT>0.3))./tempDenom;
    
    %%
    %calculate the ratio of time it spends in the wound after it has reached it for the first time
    for kk=1:numTracks
        handles.distanceNetwork.inWoundRatio(kk) = mean(handles.inWoundNetwork(1:handles.distanceNetwork.numHops(kk)-handles.distanceNetwork.numHopsTranslation(kk),kk));
    end
    %calculate the ratio of time it spends in the exploration (reaches wound and continues or not moving)
    %%    
    tempDenom                                   = handles.distanceNetwork.numHopsTranslation + handles.distanceNetwork.numHopsExploration;
    handles.distanceNetwork.inExplorationRatio  = handles.distanceNetwork.numHopsExploration./tempDenom;
    %%
    %calculate the ratio of time the neutrophil is idle once it is in the wound
    tempDenom                               = handles.distanceNetwork.numHopsExploration;

    tempNumer                               = tempDenom -sum(handles.distanceNetwork.perHopE>0.5);
    tempDenom(tempDenom==0)                 = inf;
    handles.distanceNetwork.IdleWoundRatio  = (tempNumer)./tempDenom;

    %calculate the ratio of time the neutrophil is leaving the wound 
    % LeaveWoundRatio is just relative to non-idle movements
    tempDenom                               = handles.distanceNetwork.numHopsExploration;
    tempDenom                               = tempDenom - tempNumer;
    tempDenom(tempDenom==0)                 = inf;
    tempNumer                               = sum((handles.distanceNetwork.perHopEffectiveE<-0)&(handles.distanceNetwork.perHopE>0.5));
    handles.distanceNetwork.LeaveWoundRatio = tempNumer./tempDenom;

    % LeaveWoundRatio2 is just relative to all movements
    tempDenom                               = handles.distanceNetwork.numHopsExploration;
    tempDenom(tempDenom==0)                 = inf;
    handles.distanceNetwork.LeaveWoundRatio2= tempNumer./tempDenom;

%%
else
    %If wound region is not received, then just calculate the oriented hops

    %For every track, calculate angle between initial and final points
    for kk = 1:numTracks
        ttt                                                         = 1:handles.distanceNetwork.numHops(kk);
        initFinCoord                                                = handles.nodeNetwork(handles.finalNetwork((ttt),kk),(1:3));

        tempAngle                                                   = angle((initFinCoord(end,1)-initFinCoord(1,1))+i*(initFinCoord(end,2)-initFinCoord(1,2)));
        tempXY                                                      = ([cos(tempAngle) sin(tempAngle);-sin(tempAngle) cos(tempAngle)]*initFinCoord(:,[1 2])')';
        tempD                                                       = diff(tempXY); %#ok<UDIM>

        handles.distanceNetwork.perHopOriented(ttt(1:end-1),kk)     = tempD(:,1);
        handles.distanceNetwork.perHopEffective(ttt(1:end-1),kk)    = handles.distanceNetwork.perHopOriented(ttt(1:end-1),kk)./handles.distanceNetwork.perHop(ttt(1:end-1),kk);
    end
    % Filter with a 3x gaussian to smooth jumps

    handles.distanceNetwork.perHop2             =   [handles.distanceNetwork.perHop(1,:) ;
        0.25*handles.distanceNetwork.perHop(1:end-2,:) + ...
        0.5*handles.distanceNetwork.perHop(2:end-1,:) + ...
        0.25*handles.distanceNetwork.perHop(3:end,:) ;handles.distanceNetwork.perHop(end,:)];

    handles.distanceNetwork.perHopEffective2    =   [handles.distanceNetwork.perHopEffective(1,:) ;
        0.25*handles.distanceNetwork.perHopEffective(1:end-2,:) + ...
        0.5*handles.distanceNetwork.perHopEffective(2:end-1,:) + ...
        0.25*handles.distanceNetwork.perHopEffective(3:end,:); handles.distanceNetwork.perHopEffective(end,:)];

end

%% Finally calculate the meandering ratio and average depth of the tracks if not calculated previously
    distExtBranch(1,numTracks)=0;
    averageDepthBranch(1,numTracks)=0;
    anglFromStart(1,numTracks)=0;
    
    for k=1:numTracks
        indexN                                      = handles.finalNetwork(:,k);indexN(indexN==0)=[];
        %calculate now the distance between the extremes of the branch
        distExtBranch(1,k)                          = sqrt(sum((diff(handles.nodeNetwork(indexN([1 end]),1:3))).^2,2));    
        averageDepthBranch(1,k)                     = mean(handles.nodeNetwork(indexN,3));
        
        diffBetPoints                       = diff(handles.nodeNetwork(indexN([1 end]),1:2));
        distFromStart                          = sqrt((sum(diffBetPoints.^2,2)));
        
        anglFromStart(1,k)                       = acos(diffBetPoints(:,1)./(distFromStart+1e-30));
       
    end
    % Tortuosity and meanderRatio are inverse of each other
    handles.distanceNetwork.meanderRatio            = distExtBranch./handles.distanceNetwork.totPerTrack;
    handles.distanceNetwork.tortuosity              = handles.distanceNetwork.totPerTrack./distExtBranch;

    handles.distanceNetwork.averageDepthTrack       = averageDepthBranch;
    handles.distanceNetwork.angleTrack              = anglFromStart;

