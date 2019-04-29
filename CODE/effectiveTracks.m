function handles = effectiveTracks(handles,woundRegion,thresEffecHop)
%function handles = effectiveTracks(handles,woundRegion,thresEffecHop)
%
%--------------------------------------------------------------------------
% effectiveTracks  calculates a series of maps with metrics based on the 
%     metrics previously calculated.
%
%       INPUT
%         handles:          handles struct containing distanceNetwork.
%         woundRegion:      matrix indicating location of wound area.
%         thresEffecHop:    threshold
%
%       OUTPUT
%         handles:          handles with the following new fields in 
%                           handles.distMaps:
%                           absDistPerHop     anglePerHop     oriDistPerHop   
%                           furthestPoint     furthestPoint2  latDistPerHop
%                           effDistPerHop     nodeNetwork     absDistMap
%                           oriDistMap        angleMap        latDistMap:   
%                           effDistMap
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



handles = removeHandlesFields(handles);

%%


if ~exist('thresEffecHop','var')
    thresEffecHop=0.3;
end
%first Detect the region where the wound exists
if ~exist('woundRegion','var')
        indexWound              = [];
elseif isempty(woundRegion)
        indexWound              = [];
else
    indexWound              = find (woundRegion);
end
tempIndex               = sub2ind([handles.rows handles.cols],round(handles.nodeNetwork(:,1)),round(handles.nodeNetwork(:,2)));
objectInWound           = ismember(tempIndex,indexWound);

% Get the situation of inwound(1/0) arranged in the same format as finalNetwork
q                       =handles.finalNetwork(:);

%
q2                      = find(handles.finalNetwork);
q3                      = q(q2);
q4                      = objectInWound(q3);

%
[longestTrack,numTracks]= size(handles.finalNetwork);
q(q2,2)                 = q4;
%%
handles.inWound                             = reshape(q(:,2),[longestTrack,numTracks]);     %1/0 if is in the wound (may go out)
posFinalNetwork                             = (handles.finalNetwork>0);                     %1/0 if there is an object in track
networkReachWound                           = posFinalNetwork.*(cumsum(handles.inWound)>0);                  %1/0 if the object REACHES the wound



handles.distanceNetwork.numHops             = sum(posFinalNetwork)-1;
handles.distanceNetwork.numHopsExploration  = sum(networkReachWound)-1;
handles.distanceNetwork.numHopsExploration    (handles.distanceNetwork.numHopsExploration   ==-1)=0;
handles.distanceNetwork.numHopsTranslation  = handles.distanceNetwork.numHops-handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.numHopsInWound      = sum(handles.inWound)-1;
handles.distanceNetwork.numHopsInWound        (handles.distanceNetwork.numHopsInWound   ==-1)=0;


handles                                     = trackOrientation (handles,0);
%%
avDisHop                                    = mean(abs(handles.distMaps.oriDistPerHop(handles.finalNetwork>0)));

handles.distanceNetwork.forwardRatioTot     = sum(handles.distMaps.effDistPerHop>0.6)./                                     handles.distanceNetwork.numHops;
handles.distanceNetwork.forwardRatioTot2    = sum(handles.distMaps.oriDistPerHop>0.7)./                                     handles.distanceNetwork.numHops;
handles.distanceNetwork.forwardRatioE       = sum((handles.inWound.*handles.distMaps.effDistPerHop)>0.6)./                  handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.forwardRatioT       = sum(((posFinalNetwork-networkReachWound).*handles.distMaps.effDistPerHop)>0.6)./handles.distanceNetwork.numHopsTranslation;
handles.distanceNetwork.inExplorationRatio  = handles.distanceNetwork.numHopsInWound./                                      handles.distanceNetwork.numHops;
handles.distanceNetwork.inWoundRatio        = handles.distanceNetwork.numHopsInWound./                                      handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.IdleWoundRatio      = sum(handles.inWound.*(handles.distMaps.absDistPerHop<0.7))./                  handles.distanceNetwork.numHopsInWound;
handles.distanceNetwork.backwardRatioTot    = sum(handles.distMaps.oriDistPerHop<-0.6)./                          handles.distanceNetwork.numHops;
handles.distanceNetwork.LeaveWoundRatio     = sum(networkReachWound.*(handles.distMaps.oriDistPerHop<-1)  ) ./              handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.LeaveWoundRatio2    = sum(networkReachWound.*(handles.distMaps.effDistPerHop<-0)  ) ./              handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.LeaveWoundRatio3    = sum(networkReachWound.*(handles.distMaps.oriDistPerHop<-1)  ) ./              sum(networkReachWound.*(handles.distMaps.absDistPerHop>1)  );
handles.distanceNetwork.staticRatio         = (sum(handles.distMaps.furthestPoint2)-1)                      ./              handles.distanceNetwork.numHops;
handles.distanceNetwork.returnRatio         = sum(handles.distMaps.furthestPoint(2:end,:).*(handles.distMaps.oriDistPerHop(1:end-1,:)<0)  ) ./ sum(handles.distMaps.furthestPoint(2:end,:));
handles.distanceNetwork.returnRatio2        = sum(handles.distMaps.furthestPoint(2:end,:).*(handles.distMaps.oriDistPerHop(1:end-1,:)<-(3.5*avDisHop))  ) ./ sum(handles.distMaps.furthestPoint(2:end,:));


%%
handles.distanceNetwork.absVelocity         = sum(handles.distMaps.absDistPerHop)./handles.distanceNetwork.numHops;
handles.distanceNetwork.absVelocityExp      = sum(networkReachWound.*handles.distMaps.absDistPerHop)./handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.absVelocityTra      = sum((1-networkReachWound).*handles.distMaps.absDistPerHop)./handles.distanceNetwork.numHopsTranslation;


handles.distanceNetwork.oriVelocity         = sum(handles.distMaps.oriDistPerHop)./handles.distanceNetwork.numHops;
handles.distanceNetwork.oriVelocityExp      = sum(networkReachWound.*handles.distMaps.oriDistPerHop)./handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.oriVelocityTra      = sum((1-networkReachWound).*handles.distMaps.oriDistPerHop)./handles.distanceNetwork.numHopsTranslation;


handles.distanceNetwork.latVelocity         = sum(handles.distMaps.latDistPerHop)./handles.distanceNetwork.numHops;
handles.distanceNetwork.latVelocityExp      = sum(networkReachWound.*handles.distMaps.latDistPerHop)./handles.distanceNetwork.numHopsExploration;
handles.distanceNetwork.latVelocityTra      = sum((1-networkReachWound).*handles.distMaps.latDistPerHop)./handles.distanceNetwork.numHopsTranslation;
%%
handles.distanceNetwork.crossWoundBack      = (sum(networkReachWound~=handles.inWound)>0);
handles.distanceNetwork.crossWoundBackRatio = handles.distanceNetwork.crossWoundBack/numTracks;


%% New metrics based on the activation of the neutrophils

%first find the initial and final times of the track
handles.distMaps.initialTime                = handles.nodeNetwork(handles.finalNetwork(1,:),5)';
handles.distMaps.finalTime                  = handles.distMaps.initialTime+handles.distanceNetwork.numHops;

%detect which neutrophils are active
indexActivatedNeuts                         = find(handles.distMaps.activationTime3);
handles.distMaps.numActiveNeuts             = numel(indexActivatedNeuts);

startTime                                   = handles.distMaps.activationTime3(indexActivatedNeuts) -handles.distMaps.initialTime(indexActivatedNeuts);
endTime                                     = handles.distMaps.enterWoundTime3(indexActivatedNeuts) -handles.distMaps.initialTime(indexActivatedNeuts);
endTime2                                    = handles.distMaps.finalTime(indexActivatedNeuts)       -handles.distMaps.initialTime(indexActivatedNeuts);
endTime(endTime<0)                          = endTime2(endTime<0);

%%
%startTime                                   = 
for counterActiveN = 1:handles.distMaps.numActiveNeuts
    indexTimes                              = max(1,startTime(counterActiveN)):endTime(counterActiveN);
    handles.distanceNetwork.absVelocityAC(indexActivatedNeuts(counterActiveN))     = mean(handles.distMaps.absDistPerHop(indexTimes,indexActivatedNeuts(counterActiveN)));
    handles.distanceNetwork.oriVelocityAC(indexActivatedNeuts(counterActiveN))     = mean(handles.distMaps.oriDistPerHop(indexTimes,indexActivatedNeuts(counterActiveN)));
    handles.distanceNetwork.latVelocityAC(indexActivatedNeuts(counterActiveN))     = mean(abs(handles.distMaps.latDistPerHop(indexTimes,indexActivatedNeuts(counterActiveN))));

end

% set the values of non active neutrophils to nan
indexNONActivatedNeuts                                              = find(handles.distMaps.activationTime3==0);
handles.distanceNetwork.absVelocityAC(indexNONActivatedNeuts)       = nan;
handles.distanceNetwork.oriVelocityAC(indexNONActivatedNeuts)       = nan;
handles.distanceNetwork.latVelocityAC(indexNONActivatedNeuts)       = nan;

%% get the ratio of the number of hops that are inactive vs those after activation of the neutrophil
% this is more reliable than the number of neutrophil active relative to total as small tracks would bias
% to lower ratios. To calculate the ratio of the active region only (no wound hops) change endTime2 for
% endTime

handles.distMaps.relativeActiveHops         = sum(-startTime+endTime2)/sum(handles.distanceNetwork.numHops);

%%
%[rowsData,colsData,levsData]= size(woundRegion);
rowsData    = handles.rows;
colsData    = handles.cols;
levsData    = handles.levs;

tempMaps = metricPositionMap([handles.distMaps.nodeNetwork(:,[1 2 ]) (handles.distMaps.nodeNetwork(:, 33:37))]);
[rowsMaps,colsMaps,levsMaps]= size(tempMaps);
%%
tempMaps2                           =zeros(size(tempMaps));
if handles.numFrames >200
    filtG                               =gaussF(9,9,1);
elseif handles.numFrames<100
    filtG                               =gaussF(15,15,1);
else
    filtG                               =gaussF(13,13,1);
end
    
for counterMap = 1:5
    tempMaps2(:,:,counterMap)   = imfilter(tempMaps(:,:,counterMap),filtG,'replicate');
    tempMaps2(:,:,counterMap)   = max(max(tempMaps(:,:,counterMap)))*((tempMaps2(:,:,counterMap)))/max(max(tempMaps2(:,:,counterMap))) ;
end
%%
if rowsMaps==rowsData
    handles.distMaps.absDistMap   = tempMaps2(:,:,1); 
    handles.distMaps.oriDistMap   = tempMaps2(:,:,2); 
    handles.distMaps.angleMap     = tempMaps2(:,:,3); 
    handles.distMaps.latDistMap   = tempMaps2(:,:,4); 
    handles.distMaps.effDistMap   = tempMaps2(:,:,5); 
else
    dimsR = 1:min(rowsData,rowsMaps);
    dimsC = 1:min(colsData,colsMaps);
    handles.distMaps.absDistMap   = tempMaps2(dimsR,dimsC,1); 
    handles.distMaps.oriDistMap   = tempMaps2(dimsR,dimsC,2); 
    handles.distMaps.angleMap     = tempMaps2(dimsR,dimsC,3); 
    handles.distMaps.latDistMap   = tempMaps2(dimsR,dimsC,4); 
    handles.distMaps.effDistMap   = tempMaps2(dimsR,dimsC,5);     
%     handles.distMaps.absDistMap   = tempMaps2(1:rowsData,1:colsData,1); 
%     handles.distMaps.oriDistMap   = tempMaps2(1:rowsData,1:colsData,2); 
%     handles.distMaps.angleMap     = tempMaps2(1:rowsData,1:colsData,3); 
%     handles.distMaps.latDistMap   = tempMaps2(1:rowsData,1:colsData,4); 
%     handles.distMaps.effDistMap   = tempMaps2(1:rowsData,1:colsData,5);     
end
