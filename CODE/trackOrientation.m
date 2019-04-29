function handles = trackOrientation (handles,getMaps)
%function handles = trackOrientation (handles,getMaps,numDilatation)
%
%--------------------------------------------------------------------------
% trackOrientation  calculate the track Orientation relative to all tracks 
%
%       INPUT
%         handles:      handles struct with tracks data
%         getMaps:      1 to recalculate the maps of metrics, [0] do not recalculate
%
%       OUTPUT
%         handles:      updated handles struct with all the distMaps options recalculated
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

%disable display option
if ~exist('getMaps','var')
    getMaps=0;
end


% initialise the new variables to zero by setting to zero the extreme value of each matrix
[longestTrack,numTracks]                                = size(handles.finalNetwork);

handles.distMaps.absDistPerHop (longestTrack,numTracks) = 0;
handles.distMaps.anglePerHop   (longestTrack,numTracks) = 0;
handles.distMaps.oriDistPerHop (longestTrack,numTracks) = 0;
handles.distMaps.furthestPoint (longestTrack,numTracks) = 0;
handles.distMaps.furthestPoint2 (longestTrack,numTracks) = 0;
handles.distMaps.latDistPerHop (longestTrack,numTracks) = 0;

% traverse all the tracks and calculate the absolute distance and the angle of every hop
for counterTrack = 1:numTracks
    indexN                                              = handles.finalNetwork(:,counterTrack);
    indexN(indexN==0)                                   = [];
    tempCoords                                          = handles.nodeNetwork(indexN,1:3);
    handles.distMaps.absDistPerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack)            = sqrt(sum((diff(tempCoords)).^2,2));
    handles.distMaps.anglePerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack)              = angle(diff(tempCoords(:,1)+i*tempCoords(:,2)));
end

% obtain the histogram of the angle of all those displacements larger than 0
[yy,xx]=hist(handles.distMaps.anglePerHop(handles.distMaps.absDistPerHop>0),20);

% ---- The main orientation of the MOVEMENT of the neutrophils is obtained from the ANGLE of all tracks --- 
% the mean is required in case the maximum is not a single point
MainOrientation = mean(xx(yy==max(yy)));
% a restriction into the positive side of X for cases where the cells are moving towards the right


% The coordinates of the system are rotated so that the oriented distance is in one axis and the lateral
% on the opposite and then the calculation is simplified

RotatedCoordinates                                     = ([cos(MainOrientation) sin(MainOrientation);-sin(MainOrientation) cos(MainOrientation)]*handles.nodeNetwork(:,[1 2])')';

ExtremeRotatedCoordinates                               = max(RotatedCoordinates);

%%
avDisHop                                                            = mean(handles.distMaps.absDistPerHop(handles.finalNetwork>0));
tempCoords3(handles.numFrames,1)                                    = 0;
tempCoords4(handles.numFrames,1)                                    = 0;

%%          
for counterTrack = 1:numTracks
    % again traverse the tracks to extract the parameters from each from the new coordinates
    indexN                                                          = handles.finalNetwork(:,counterTrack);
    indexN(indexN==0)                                               = [];

    lengthCurrTrack                                                 = size(indexN,1);
    tempCoords                                                      = RotatedCoordinates(indexN,1:2);
    tempCoords2                                                     = handles.nodeNetwork(indexN,1:2);


    [maxDispl,indexMaxDisp]                                         = max(tempCoords(:,1));
    [maxDispl2,indexMaxDisp2]                                       = max(tempCoords2(:,2));

    finalCurrentRow                                                 = handles.distanceNetwork.numHops(counterTrack);
    handles.distMaps.oriDistPerHop(1:finalCurrentRow,counterTrack)            = diff(tempCoords(:,1));
    handles.distMaps.latDistPerHop(1:finalCurrentRow,counterTrack)            = diff(tempCoords(:,2));
    handles.distMaps.furthestPoint(indexMaxDisp:1+finalCurrentRow,counterTrack) = 1;
    handles.distMaps.furthestPoint2(tempCoords(:,1)>(maxDispl-(8*avDisHop)),counterTrack) = 1;
    handles.distMaps.furthestPoint3(1:2,counterTrack)               = [maxDispl2;tempCoords2(indexMaxDisp2,1)];
    handles.distMaps.initialPoint3(1:2,counterTrack)                = [tempCoords2(1,2);tempCoords2(1,1)];
    handles.distMaps.initialPoint2(1:2,counterTrack)                = [tempCoords(1,2);tempCoords(1,1)];
    %% tempCoords3 is the oriented distance of the current track with a 5x low pass filter


    currOrientedDistance                                            =  handles.distMaps.oriDistPerHop(:,counterTrack);
    tempCoords3                                                     =   0.14*currOrientedDistance(1:end-4)+...
                                                                        0.23*currOrientedDistance(2:end-3)+...
                                                                        0.26*currOrientedDistance(3:end-2)+...
                                                                        0.23*currOrientedDistance(4:end-1)+...
                                                                        0.14*currOrientedDistance(5:end);   
                                                                
    tempCoords3                                                     = [tempCoords3(1);tempCoords3(1);tempCoords3;tempCoords3(end);tempCoords3(end)] ;    
    tempCoords3(handles.numFrames,1)                                = 0;
    
    %find if the neutrophil gets to the wound or not
    isCurrNeutInWound                                               = find(handles.inWound(:,counterTrack)>0,1);
    %% tempCoords4 looks for the next (40 frames OR point where track ends OR enters wound) and gets the mean of highest 80% hops 
    %clear tempCoords4
    for counterPosition =1:handles.numFrames
        lastPointCurrData                                           = min([handles.numFrames counterPosition+40 lengthCurrTrack isCurrNeutInWound+5]);
        tempCoords5                                                 = sort(currOrientedDistance(counterPosition:lastPointCurrData));
        if ~isempty(tempCoords5)
            oneFiftheTempCoords                                     = ceil(numel(tempCoords5)/5);
            tempCoords4(counterPosition,1)                          = mean(tempCoords5(oneFiftheTempCoords:end));
        else
            tempCoords4(counterPosition,1)                          = 0;
        end

    end
                                                                
    %%
    %assign the point where the track enters the wound
    if isempty(isCurrNeutInWound)
        handles.distMaps.enterWoundPoint3(1:2,counterTrack)         = [0;0];
        handles.distMaps.enterWoundTime3(1,counterTrack)            = 0;
        
    else
        handles.distMaps.enterWoundPoint3(1:2,counterTrack)         = [tempCoords2(isCurrNeutInWound,2);tempCoords2(isCurrNeutInWound,1)];
        handles.distMaps.enterWoundTime3(1,counterTrack)            = isCurrNeutInWound+handles.nodeNetwork(handles.finalNetwork(1,counterTrack),5);
        
    end
    % detect if the neutrophil starts movement towards the wound and when it starts
    isCurrNeutActive                                                = find((tempCoords3>1)&(tempCoords4>1),1);

    if isempty(isCurrNeutActive)
        handles.distMaps.activationPoint3(1:2,counterTrack)         = [0;0];
        handles.distMaps.activationTime3(1,counterTrack)            = 0;
        handles.distMaps.activationDistTW(1,counterTrack)           = 0;

    else
        %to be an active neutrophil, the track has to be AT LEAST, 15 hops in length after the initiation
        if ((lengthCurrTrack-isCurrNeutActive)>15)
            if (isCurrNeutActive<isCurrNeutInWound)
                handles.distMaps.activationPoint3(1:2,counterTrack)     = [tempCoords2(isCurrNeutActive,2);tempCoords2(isCurrNeutActive,1)];
                handles.distMaps.activationTime3(1,counterTrack)        = isCurrNeutActive+handles.nodeNetwork(handles.finalNetwork(1,counterTrack),5);
                handles.distMaps.activationDistTW(1,counterTrack)       = ExtremeRotatedCoordinates(1) - handles.distMaps.initialPoint2(2,counterTrack);
            else
                %if it DOES reach the wound BEFORE initiation, speed high for at least 5 hops
                if (mean(tempCoords3(isCurrNeutActive:isCurrNeutActive+5))>1)
                    handles.distMaps.activationPoint3(1:2,counterTrack)     = [tempCoords2(isCurrNeutActive,2);tempCoords2(isCurrNeutActive,1)];
                    handles.distMaps.activationTime3(1,counterTrack)        = isCurrNeutActive+handles.nodeNetwork(handles.finalNetwork(1,counterTrack),5);
                    handles.distMaps.activationDistTW(1,counterTrack)       = ExtremeRotatedCoordinates(1) - handles.distMaps.initialPoint2(2,counterTrack);


                else
                    handles.distMaps.activationPoint3(1:2,counterTrack)         = [0;0];
                    handles.distMaps.activationTime3(1,counterTrack)            = 0;
                    handles.distMaps.activationDistTW(1,counterTrack)           = 0;


                end
            end        
        else
            handles.distMaps.activationPoint3(1:2,counterTrack)         = [0;0];
            handles.distMaps.activationTime3(1,counterTrack)            = 0;
            handles.distMaps.activationDistTW(1,counterTrack)           = 0;

        end
    end
    
end
%%  

handles.distMaps.effDistPerHop                          = handles.distMaps.oriDistPerHop./handles.distMaps.absDistPerHop;
handles.distMaps.effDistPerHop(isnan(handles.distMaps.effDistPerHop)) = 0;
handles.distMaps.nodeNetwork(:,[1:3 4 9 10 22:27]) = handles.nodeNetwork(:,[1:3 4 9 10 22:27]);


for counterTrack = 1:numTracks
%
    indexN                  =handles.finalNetwork(:,counterTrack);
    indexN(indexN==0)       =[];

    handles.distMaps.nodeNetwork(indexN(2:end),33) =     handles.distMaps.absDistPerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack);
    handles.distMaps.nodeNetwork(indexN(2:end),34) =     handles.distMaps.oriDistPerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack);
    handles.distMaps.nodeNetwork(indexN(2:end),35) =     handles.distMaps.anglePerHop (1:handles.distanceNetwork.numHops(counterTrack),counterTrack);
    handles.distMaps.nodeNetwork(indexN(2:end),36) = abs(handles.distMaps.latDistPerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack));
    handles.distMaps.nodeNetwork(indexN(2:end),37) =     handles.distMaps.effDistPerHop(1:handles.distanceNetwork.numHops(counterTrack),counterTrack);
%
end
%% Get Maps of metrics
if getMaps==1

    tempMaps = metricPositionMap([handles.distMaps.nodeNetwork(:,[1 2 ]) (handles.distMaps.nodeNetwork(:, 33:37))]);
    handles.distMaps.absDistMap   = tempMaps(:,:,1); 
    handles.distMaps.oriDistMap   = tempMaps(:,:,2); 
    handles.distMaps.angleMap     = tempMaps(:,:,3); 
    handles.distMaps.latDistMap   = tempMaps(:,:,4); 
    handles.distMaps.effDistMap   = tempMaps(:,:,5); 
end

