function [hPlotNet,distBetweenTracks] =plotTrackStats(handles,typeOfPlot,currentTrack,micronsPerPixel,framesPerSec,secondTrack)
%function hPlotNet=plotTrackStats(handles)
%function hPlotNet=plotTrackStats(handles,typeOfPlot)
%function
%hPlotNet=plotTrackStats(handles,typeOfPlot,currentTrack,...
%                        micronsPerPixel,framesPerSec)
%
%--------------------------------------------------------------------------
% plotTrackStats  Advance command to plot the STATS from the tracks of a
%     neutrophil analysis, it plots in many ways
%
%       INPUT
%         handles: handles struct containing:
%                        handles.nodeNetwork   [numRBC detected x 12 params]
%                      handles. finalNetwork   either 1 track or
%                      [depth of tracks x numTracks]
%
%         typeOfPlot:           a number with the following options:
%                                       0  = distance per frame      single track
%                                       1  = angle per frame         single track
%                                       2  = distance from start     single track
%                                       3  = angle from start        single track
%                                       4  = cumulative distance     single track
%                                       5  = distance per frame      all tracks
%                                       6  = angle per frame         all tracks
%                                       7  = distance from start     all tracks
%                                       8  = angle from start        all tracks
%                                       9  = cumulative distance     all tracks
%                                       10 = av distance vs angle rose plot
%                                       11 = av velocity vs angle rose plot
%                                       11 = Total track distance as a polar plot
%                                       12 = Maximum track velocity as a polar plot
%                                       13 = meandering index (ratio) as a polar plot
%                                       14 = the tracks in velocity
%                                       and time as "events"
%
%         microsPerPixel:       the calculations of distance and
%         velocity are made on a
%
%         framesPerPixel:       pixel/frame basis unless there is a
%         distanceProp to give proportionallity
%
%         currentTrack:         a single track to plot, otherwise
%         the first one is selected
%
%       OUTPUT
%         hPlotNet:             a handle to the plot with all the results
%
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
%------ no input data is received, error -------------------------
%------ at least 2 parameters are required
if nargin <1 ;
    help plotTrackStats;
    hPlotNet        = [];
    return;
else
    nodeNetwork     = handles.nodeNetwork(:,[1 2 5]);
    finalNetwork    = handles.finalNetwork;
end
numTracks           = size(finalNetwork,2);
%%
%------ arguments received revision   ----------------------------
if ~exist('typeOfPlot','var')
    typeOfPlot      =1;
end
if (typeOfPlot<-1)||(typeOfPlot>15)
    help plotTrackStats;
    hPlotNet        =[];
    return;
end

if ~exist('currentTrack','var')
    currentTrack=1;
end
if (exist('micronsPerPixel','var'))&(isempty(micronsPerPixel))
    clear micronsPerPixel;
end
if (exist('framesPerSec','var'))&(isempty(framesPerSec))
    clear framesPerSec;
end


%define all the labels according to the input arguments
if (~exist('micronsPerPixel','var'))||(~exist('framesPerSec','var'))
    lab_Velo_mm_s  = 'Velocity [pix/frame]';
else
    lab_Velo_mm_s  = 'Velocity [\mum /s]';
end

if ~exist('micronsPerPixel','var') ;
    micronsPerPixel = 1;
    lab_Cols_dist = 'Columns [pixels]';
    lab_Rows_dist = 'Rows [pixels]';
    yLabUnits =' [pixels]';
else
    lab_Rows_dist = 'Rows [ \mum ]';
    lab_Cols_dist = 'Columns [ \mum ]';
    yLabUnits = ' [\mu m]';
end

if ~exist('framesPerSec','var') ;
    framesPerSec = 1;
    lab_Time_fps = 'Time [frames]'  ;
else
    lab_Time_fps = 'Time [sec]';
end
%%
%------ regular size check ---------------------------------------
[numNodes,numCols]          = size(handles.nodeNetwork);
[maxDepthTrack,numTracks]   = size(handles.finalNetwork);
%%
% Get axis ready
h0                          = gcf;
h1                          = gca;
hChild                      = allchild(h0);
hColP                       = findobj(hChild,'tag','colorBarPlot');
hColB                       = findobj(hChild,'tag','colorbarLabel');
set(hColP,'visible','off');
set(hColB,'visible','off');
hold on;
set(h1,'Color',[1 1 1]);


%------- colour codes for the plot
colorID2 = {'o', 'x', 's', '+', 'v',  '*', '^', '<', '>','h','p','o'};
colorID3=0.9*[ [0 0 1];[0 1 0];[1 0 0];...
    [0 0 0.75]; [0 0 0.25]; [0 0 0];...
    [0 1 0];    [0 0.75 0]; [0 0.5 0];...
    [1 0 0];    [0.75 0 0]; [0.35 0 0];...
    [0 1 1];    [1 0 1];    [1 1 0 ];...
    [0.5 0.5 0];[0.5 0 0.5];[0 0.5 0.5];...
    [0.25 0.25 0.25];[0.25 1 0];[0.25 0 0.25];[0 0.25 0.25]];

colorID6                            = interp1(jet,linspace(1,64,numTracks));

if ~isfield(handles,'meanderingRatio')
    handles=effectiveDistance(handles);
end
%-------- The structure of nodeNetwork will be [X Y z] each column
%-------- is separately transformed int
%-------- a matrix [depthB x numTracks] which will be used as X,Y,Z
%for the plot3 function
%-------- how to plot some nets that DO NOT have same number of nodes????
nodeNetworkX = (handles.nodeNetwork(:,2));
nodeNetworkY = (handles.nodeNetwork(:,1));
nodeNetworkZ = (handles.nodeNetwork(:,5));

%----- a neighbourhood for the group of tracks has been received
neighNetwork = handles.finalNetwork;
%----- finalNetwork will keep a vector with the lengths of the branches
hPlotNet.finalNetwork = sum(neighNetwork>0);
hPlotNet.numTracks = numTracks;
hPlotNet.relDistPrevHop = [];
hPlotNet.relAnglPrevHop = [];
 
if typeOfPlot<5
    % all plots here are for a single track, the plotting points are selected
    % according to currentTrack and then XX,YY,ZZ are obtained
    if typeOfPlot==-1
        if (currentTrack>numTracks)|(secondTrack>numTracks)
            disp('The tracks selected are not valid.')
            hPlotNet        = [];distBetweenTracks=[];
            return;
        end
        
        plottingPoints_1    = neighNetwork(1:hPlotNet.finalNetwork(currentTrack),currentTrack);
        plottingPoints_2    = neighNetwork(1:hPlotNet.finalNetwork(currentTrack),secondTrack);
        if numel(plottingPoints_1)<3
            disp('Track has less than 3 points, no plots will be generated')
            hPlotNet        = [];distBetweenTracks=[];
            return;
        end
        % Find the valid positions ON THE neighNetwork
        posTrack_1          = find(plottingPoints_1);
        posTrack_2          = find(plottingPoints_2);
        
        ZZ_1                = nodeNetworkZ(plottingPoints_1(posTrack_1));
        ZZ_2                = nodeNetworkZ(plottingPoints_2(posTrack_2));

        posOverlap          = intersect(ZZ_1,ZZ_2);
        posTrack_1B         = ismember(ZZ_1,posOverlap);
        posTrack_2B         = ismember(ZZ_2,posOverlap);
        
        XX_1                = nodeNetworkX(plottingPoints_1(posTrack_1B));
        YY_1                = nodeNetworkY(plottingPoints_1(posTrack_1B));

        XX_2                = nodeNetworkX(plottingPoints_2(posTrack_2B));
        YY_2                = nodeNetworkY(plottingPoints_2(posTrack_2B));
        
        distBetweenTracks(1:handles.numFrames)   = -1;
        distBetweenTracks(posOverlap)   = sqrt((XX_1-XX_2).^2 +(YY_1-YY_2).^2  );
        
        jumpsTime  = find(diff(posOverlap)>1);
        
        if (any(jumpsTime))
            %these positions will be empty
            emptyPositions = posOverlap(jumpsTime)+1;
            distBetweenTracks(emptyPositions) = 0.5*distBetweenTracks(emptyPositions-1)+0.5*distBetweenTracks(emptyPositions+1);
        end
        
        %-----distance between two tracks
     
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(1:handles.numFrames,distBetweenTracks,'marker',...
                colorID2{1+rem(currentTrack-1,12)},...
                'color',colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis tight;grid on;axis xy;
            ylabel(strcat('Distance between Tracks',yLabUnits),'fontsize',11);

        
        
        
    else
        plottingPoints = neighNetwork(1:hPlotNet.finalNetwork(currentTrack),currentTrack);
        if numel(plottingPoints)<3
            disp('Track has less than 3 points, no plots will be generated')
            return;
        end
        XX              = nodeNetworkX(plottingPoints);
        YY              = nodeNetworkY(plottingPoints);
        ZZ              = nodeNetworkZ(plottingPoints);
        %%
        diffBetPoints   = diff([XX YY]);
        diffBetStart    = [XX(1)-XX(2:end) YY(1)-YY(2:end)];
        distPerHop      = sqrt((sum(diffBetPoints.^2,2)));
        %anglPerHop = acos(diffBetPoints(:,1)./(distPerHop+1e-30));
        distFromStart   = sqrt(sum(diffBetStart.^2,2));
        anglFromStart   = acos(diffBetStart(:,1)./(distFromStart+1e-30));
        
        %The angles require a complicated calculation to get it right
        diffBet2Points2 = (-[XX(1:end-2) YY(1:end-2)]+[XX(3:end) YY(3:end)]);
        distPer2Hop     = sqrt((sum(diffBet2Points2.^2,2)));
        distPerHop2     = (distPerHop.^2);
        distPer2Hop2    = (distPer2Hop.^2);
        distPerHop2_consec = (distPerHop2(1:end-1)+distPerHop2(2:end));
        if ~isempty(distPer2Hop2)
            q=[0;0;distPerHop2_consec==distPer2Hop2];
        else
            q=0;
        end
        if any(q)
            %q2=(distPerHop2_consec==distPer2Hop2);
            XX(q==1) = XX(q==1)+1.234e-3;
            YY(q==1) = YY(q==1)+1.432e-3;
            %distPerHop2_consec(q2==1) = distPerHop2_consec(q2==1)+1.432e-3;
            diffBetPoints = diff([XX YY]);
            distPerHop = sqrt((sum(diffBetPoints.^2,2)));
            diffBet2Points2 = (-[XX(1:end-2) YY(1:end-2)]+[XX(3:end) YY(3:end)]);
            distPer2Hop = sqrt((sum(diffBet2Points2.^2,2)));
            distPerHop2 = (distPerHop.^2);
            distPer2Hop2 = (distPer2Hop.^2);
            distPerHop2_consec = (distPerHop2(1:end-1)+distPerHop2(2:end));
        end
        for counterP=1:size(YY,1)-2
            theta_a = -atan((YY(counterP+1)-YY(counterP))/(XX(counterP+1)-XX(counterP)));
            rotMat =  [cos(theta_a) -sin(theta_a);sin(theta_a) cos(theta_a)];
            p2 = rotMat*[XX(counterP+1);YY(counterP+1)];
            p3 = rotMat*[XX(counterP+2);YY(counterP+2)];
            p33 = p3-p2;
            anglPerHop3 = atan(p33(2)/p33(1)) ;
            signAngl3 = sign(anglPerHop3);
            if distPerHop2_consec(counterP)>distPer2Hop2(counterP)
                anglPerHop2(counterP) = -signAngl3*pi+ anglPerHop3;
            else
                anglPerHop2(counterP) = anglPerHop3;
            end
        end
        anglPerHop2(isnan(anglPerHop2)) = 0;
        anglPerHopCum       = abs(cumsum(anglPerHop2));
        anglPerHopCum    = imfilter(anglPerHopCum,[0.25 0.5 0.25],'replicate');
        %find number of turns
        currTurns = 0;
        pointsTurn =[];
        %currPoint = 1;
        for counterPoint=1:numel(anglPerHopCum)
            if (anglPerHopCum(counterPoint))>2*pi
                pointsTurn=[pointsTurn counterPoint];
                currTurns = currTurns +1;
                anglPerHopCum(counterPoint+1:end) = ...
                    anglPerHopCum(counterPoint+1:end) -2*pi;
            elseif (anglPerHopCum(counterPoint))<(-2*pi)
                pointsTurn=[pointsTurn -counterPoint];
                %disp(-counterPoint)
                currTurns = currTurns +1;
                anglPerHopCum(counterPoint+1:end) = ...
                    anglPerHopCum(counterPoint+1:end) +2*pi;
            end
        end
        %%
        distPerHop = distPerHop * micronsPerPixel;
        distFromStart = distFromStart* micronsPerPixel;
        
        if typeOfPlot==0
            %-----distance per frame
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(ZZ(2:end),distPerHop,'marker',...
                colorID2{1+rem(currentTrack-1,12)},...
                'color',colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis tight;grid on;rotate3d on;axis xy;
            ylabel(strcat('Distance per frame',yLabUnits),'fontsize',11);
        elseif typeOfPlot==1
            %-----angle per frame
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(ZZ(2:end-1),anglPerHop2,'marker',...
                colorID2{1+rem(currentTrack-1,12)},...
                'color',colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis tight;grid on;rotate3d on;axis xy;
            
            if currTurns>0
                hold on
                plot(ZZ(1)+abs(pointsTurn),sign(pointsTurn),'r*','markersize',9);
                plot(ZZ(1)+abs(pointsTurn),sign(pointsTurn),'ko','markersize',9);
                title(strcat('Number of turns = ',num2str(currTurns)),'fontsize',12);
            end
            ylabel('Angle  per frame [rad]','fontsize',11);
            
        elseif typeOfPlot==2
            %-----distance from start
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(ZZ(2:end),distFromStart,'marker',...
                colorID2{1+rem(currentTrack-1,12)},...
                'color',colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis tight;grid on;rotate3d on;axis xy;
            ylabel(strcat('Distance from origin ',yLabUnits),'fontsize',11);
        elseif typeOfPlot==3
            %-----angle from start
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(ZZ(2:end),anglFromStart,'marker',...
                colorID2{1+rem(currentTrack-1,12)},...
                'color',colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis([ZZ(1) ZZ(end) 0 pi]);grid on;rotate3d on;axis xy;
            ylabel('Angle from starting frame [rad]','fontsize',11);
        elseif typeOfPlot==4
            %-----cumulative distance traversed
            hold off;
            hPlotNet.handlePlot(currentTrack) = ...
                plot(ZZ(2:end),cumsum(distPerHop),'marker',...
                colorID2{1+rem(currentTrack-1,12)},'color',...
                colorID3(1+rem(currentTrack-1,20),:),'markersize',5,...
                'linewidth',1.5);
            axis tight;grid on;rotate3d on;axis xy;
            ylabel(strcat('Cumulative distance ',yLabUnits),'fontsize',11);
        end
    end
    xlabel(lab_Time_fps,'fontsize',18);
    xTicksTime=num2str((get(h1,'xtick')/framesPerSec)',3);
    set(h1,'xticklabel',xTicksTime);
    
elseif (typeOfPlot<10)
    % plots here are for all tracks therefore there is a loop and each track is
    % calculated and displayed
    [m1,m2] = sort(handles.distanceNetwork.totPerTrack);
    numTurnsPerTrack=[];
    
    for currentTrack = 1:numTracks
        plottingPoints = ...
            neighNetwork(1:hPlotNet.finalNetwork(currentTrack),currentTrack);
        if numel(plottingPoints)>2
            XX = nodeNetworkX(plottingPoints);
            YY = nodeNetworkY(plottingPoints);
            ZZ = nodeNetworkZ(plottingPoints);
            %%
            diffBetPoints = diff([XX YY]);
            diffBetStart = [XX(1)-XX(2:end) YY(1)-YY(2:end)];
            distPerHop = sqrt((sum(diffBetPoints.^2,2)));
            anglPerHop = acos(diffBetPoints(:,1)./(distPerHop+1e-30));
            distFromStart = sqrt(sum(diffBetStart.^2,2));
            anglFromStart = acos(diffBetStart(:,1)./(distFromStart+1e-30));
            
            distPerHop = distPerHop * micronsPerPixel;
            distFromStart = distFromStart* micronsPerPixel;
            
            %yCoord = m2(currentTrack)*ones(size(distPerHop));
            yCoord = currentTrack*ones(size(distPerHop));
            
            %---------------------------------------------------
            %The angles require a complicated calculation to get it right
            diffBet2Points2 = (-[XX(1:end-2) YY(1:end-2)]+[XX(3:end) YY(3:end)]);
            distPer2Hop = sqrt((sum(diffBet2Points2.^2,2)));
            distPerHop2 = (distPerHop.^2);
            distPer2Hop2 = (distPer2Hop.^2);
            distPerHop2_consec = (distPerHop2(1:end-1)+distPerHop2(2:end));
            if ~isempty(distPer2Hop2)
                q=[0;0;distPerHop2_consec==distPer2Hop2];
            else
                q=0;
            end
            
            %q=[0;0;distPerHop2_consec==distPer2Hop2];
            if any(q)
                %q2=(distPerHop2_consec==distPer2Hop2);
                XX(q==1) = XX(q==1)+1.234e-3;
                YY(q==1) = YY(q==1)+1.432e-3;
                %distPerHop2_consec(q2==1) = distPerHop2_consec(q2==1)+1.432e-3;
                diffBetPoints = diff([XX YY]);
                distPerHop = sqrt((sum(diffBetPoints.^2,2)));
                diffBet2Points2 = ...
                    (-[XX(1:end-2) YY(1:end-2)]+[XX(3:end) YY(3:end)]);
                distPer2Hop = sqrt((sum(diffBet2Points2.^2,2)));
                distPerHop2 = (distPerHop.^2);
                distPer2Hop2 = (distPer2Hop.^2);
                distPerHop2_consec = (distPerHop2(1:end-1)+distPerHop2(2:end));
            end
            clear anglPerHop2
            for counterP=1:size(YY,1)-2
                theta_a = -atan((YY(counterP+1)-YY(counterP))/...
                    (XX(counterP+1)-XX(counterP)));
                rotMat =  [cos(theta_a) -sin(theta_a);sin(theta_a) cos(theta_a)];
                p2 = rotMat*[XX(counterP+1);YY(counterP+1)];
                p3 = rotMat*[XX(counterP+2);YY(counterP+2)];
                p33 = p3-p2;
                anglPerHop3 = atan(p33(2)/p33(1)) ;
                signAngl3 = sign(anglPerHop3);
                if distPerHop2_consec(counterP)>distPer2Hop2(counterP)
                    anglPerHop2(counterP) = -signAngl3*pi+ anglPerHop3;
                else
                    anglPerHop2(counterP) = anglPerHop3;
                end
            end
            anglPerHop2(isnan(anglPerHop2)) = 0;
            anglPerHopCum = abs(cumsum(anglPerHop2));
            %find number of turns
            currTurns = 0;
            pointsTurn =[];
            %currPoint = 1;
            for counterPoint=1:numel(anglPerHopCum)
                if (anglPerHopCum(counterPoint))>2*pi
                    pointsTurn = [pointsTurn counterPoint];
                    currTurns = currTurns +1;
                    anglPerHopCum(counterPoint+1:end) = ...
                        anglPerHopCum(counterPoint+1:end) -2*pi;
                elseif (anglPerHopCum(counterPoint))<(-2*pi)
                    pointsTurn=[pointsTurn -counterPoint];
                    %disp(-counterPoint)
                    currTurns = currTurns +1;
                    anglPerHopCum(counterPoint+1:end) = ...
                        anglPerHopCum(counterPoint+1:end) +2*pi;
                end
                
            end
            %---------------------------------------------------
            numTurnsPerTrack(currentTrack) = currTurns;
            
            if typeOfPlot==5
                %-----distance per frame
                tempCoordValue = distPerHop;
            elseif typeOfPlot==6
                %-----angle per frame
                tempCoordValue = [0 anglPerHop2]';
                
                if currTurns>0
                    %hold on
                    plot3(ZZ(1)+abs(pointsTurn),yCoord(1),...
                        3.5,'r*','markersize',9);
                    plot3(ZZ(1)+abs(pointsTurn),yCoord(1),...
                        3.5,'ko','markersize',9);
                    %title(strcat('Number of turns = ...
                    %',num2str(currTurns)),...
                    %'fontsize',12);
                end
                
            elseif typeOfPlot==7
                %-----distance from start
                tempCoordValue = distFromStart;
            elseif typeOfPlot==8
                %-----angle from start
                tempCoordValue = anglFromStart;
            elseif typeOfPlot==9
                %-----cumulative distance traversed
                tempCoordValue = cumsum(distPerHop);
            end
            hPlotNet.handlePlot(currentTrack) = ...
                plot3(ZZ(2:end),yCoord,tempCoordValue,'marker',...
                colorID2{1+rem(currentTrack-1,11)},'color',...
                colorID6((currentTrack),:),'markersize',5,'linewidth',1.5);
        else
            numTurnsPerTrack(currentTrack) = 0;
        end
    end
    axis tight;grid on;rotate3d on;axis xy;
    xlabel(lab_Time_fps,'fontsize',14);
    ylabel('Number of Track','fontsize',14);
    xTicksTime=num2str((get(h1,'xtick')/framesPerSec)',3);
    set(h1,'xticklabel',xTicksTime);
    
    
    if typeOfPlot==5
        %-----distance per frame
        zlabel(strcat('Distance per frame',yLabUnits),'fontsize',11);
        
    elseif typeOfPlot==6
        %-----angle per frame
        ylabel('Number of Track (Num. Turns)','fontsize',14);
        zlabel('Angle  per frame [rad]','fontsize',11);
        ylabelTurns={};
        for counterTracks =1:numTracks
            ylabelTurns{counterTracks}=...
                strcat(num2str(counterTracks),'  ( ',...
                num2str(numTurnsPerTrack(counterTracks)),')');
        end
        if numTracks<15
            set(gca,'yticklabel',ylabelTurns);
            set(gca,'ytick',(1:numTracks));
        else
            tracksWithTurns=find(numTurnsPerTrack);
            if isempty(tracksWithTurns)
                tracksWithTurns =  get(gca,'ytick');
            end
            set(gca,'yticklabel',ylabelTurns(tracksWithTurns));
            set(gca,'ytick',tracksWithTurns);
        end
    elseif typeOfPlot==7
        %-----distance from start
        zlabel(strcat('Distance from origin ',yLabUnits),'fontsize',11);
    elseif typeOfPlot==8
        %-----angle from start
        zlabel('Angle from starting frame [rad]','fontsize',11);
    elseif typeOfPlot==9
        %-----cumulative distance traversed
        zlabel(strcat('Cumulative distance ',yLabUnits),'fontsize',11);
    end
    
    
    view(45,45)
elseif (typeOfPlot<14) %Eliminate unnecesary condition (typeOfPlot>9)
    %polar plots
    hold off
    if typeOfPlot==10
        %-----average distance per angle
        polar(handles.distanceNetwork.angleTrack,...
            handles.distanceNetwork.avPerTrack,'o')
        xlabel(lab_Velo_mm_s,'fontsize',14);
        ylabel('Average Track Velocity','fontsize',14);
        
    elseif typeOfPlot==11
        %--- total distance of track
        polar(handles.distanceNetwork.angleTrack,...
            handles.distanceNetwork.totPerTrack,'o')
        xlabel(lab_Velo_mm_s,'fontsize',14);
        ylabel('Total Track Distance','fontsize',14);
        
    elseif typeOfPlot==12
        %--- maximum distance of one hop
        polar(handles.distanceNetwork.angleTrack,...
            handles.distanceNetwork.maxPerTrack,'o')
        xlabel(lab_Velo_mm_s,'fontsize',14);
        ylabel('Maximum Track Velocity (largest hop)','fontsize',14);
    elseif typeOfPlot==13
        %--- meandering ratio
        polar(handles.distanceNetwork.angleTrack,...
            handles.distanceNetwork.meanderRatio,'o')
        ylabel('Meandering Ratio','fontsize',14);
    end
elseif (typeOfPlot==14)
    if numel(currentTrack) ~= numTracks
        selectSomeTracks = ones(numTracks,1);
    else
        selectSomeTracks = currentTrack;
    end
    
    for counterTrack=1:numTracks
        if selectSomeTracks(counterTrack)==1
            %---- determine WHICH nodes (points) correspond to the current Track
            plottingPoints = ...
                neighNetwork(1:hPlotNet.finalNetwork(counterTrack),...
                counterTrack);
            XX = nodeNetworkX(plottingPoints);
            YY = nodeNetworkY(plottingPoints);
            ZZ = nodeNetworkZ(plottingPoints);
            diffBetPoints = diff([XX YY]);
            diffBetStart = [XX(1)-XX(2:end) YY(1)-YY(2:end)];
            distPerHop = sqrt((sum(diffBetPoints.^2,2)));
            %distPerHop_smooth  = conv(distPerHop,gaussF(9,1,1),'same');
            anglPerHop = acos(diffBetPoints(:,1)./(distPerHop+1e-30));
            distFromStart = sqrt(sum(diffBetStart.^2,2));
            
            hPlotNet.anglFromStart(counterTrack) = ...
                atan2( (YY(end)-YY(1)) , ( (XX(end)-XX(1))+1e-30));
            hPlotNet.avDistPerTrack(counterTrack) = mean(distPerHop);
            hPlotNet.totDistPerTrack(counterTrack) = sum(distPerHop);
            hPlotNet.maxDistPerTrack(counterTrack) = max(distPerHop);
            hPlotNet.startFrame(counterTrack) = ZZ(1);
            hPlotNet.stopFrame(counterTrack) = ZZ(end);
            hPlotNet.avX(counterTrack) = mean(XX);
            hPlotNet.avY(counterTrack) = mean(YY);
            hPlotNet.anglPerHop(counterTrack) = mean(anglPerHop);
            distV=distPerHop;
            anglV=anglPerHop;
            relDistP=distV(2:end)./(distV(1:end-1)+1e-30);
            relAnglP=anglV(2:end)-anglV(1:end-1);
            relAnglP(relDistP>1000)=[];
            relDistP(relDistP>1000)=[];
            hPlotNet.relDistPrevHop = [hPlotNet.relDistPrevHop;relDistP];
            hPlotNet.relAnglPrevHop = [hPlotNet.relAnglPrevHop;relAnglP];
            
            hPlotNet.tempVelocity(counterTrack,round([framesPerSec*...
                hPlotNet.startFrame(end):framesPerSec*...
                hPlotNet.stopFrame(end)]))=...
                framesPerSec*hPlotNet.avDistPerTrack(end);
            line([hPlotNet.startFrame(end) hPlotNet.stopFrame(end)],...
                framesPerSec*hPlotNet.avDistPerTrack(end)*[1 1],...
                'marker','s','linewidth',1);
        end
    end
    
    tempAxis = [framesPerSec*min(hPlotNet.startFrame):...
        framesPerSec*max(hPlotNet.stopFrame)];
    hPlotNet.tempVelocity(isnan(hPlotNet.tempVelocity))=0;
    tempVelocityPerFrame = sum(hPlotNet.tempVelocity);
    tempEventPerFrame = sum(hPlotNet.tempVelocity>0);
    tempEventPerFrame(tempEventPerFrame==0) =   inf;
    numEventsNoZero = sum(tempEventPerFrame==inf)+1;
    temporalVelocity =  tempVelocityPerFrame./tempEventPerFrame;
    q1=cumsum(temporalVelocity);
    breakPoints=find(tempEventPerFrame==inf);
    if numel(breakPoints)>0
        temporalVelocity2=temporalVelocity;
        temporalVelocity2(breakPoints)=...
            -[q1(breakPoints(1)) diff(q1(breakPoints))];
        q3=cumsum(temporalVelocity2);
        breakPoints2=[breakPoints-1 size(q3,2)];
        breakPoints3=[breakPoints2(1) -1+diff(breakPoints2)];
        
        hPlotNet.avVelPerEvent=q3(breakPoints2)./(breakPoints3+1e-100);
        timeSlots=[1 breakPoints+1;breakPoints2];
        
        timeSlots(:,(hPlotNet.avVelPerEvent==0)|...
            (hPlotNet.avVelPerEvent>1e87)|...
            (hPlotNet.avVelPerEvent<-1e87))=[];
        
        hPlotNet.avVelPerEvent(hPlotNet.avVelPerEvent==0)=[];
        hPlotNet.avVelPerEvent((hPlotNet.avVelPerEvent>1e87))=[];
        hPlotNet.avVelPerEvent((hPlotNet.avVelPerEvent<-1e87))=[];
        
        plot(timeSlots/framesPerSec,(hPlotNet.avVelPerEvent'*[1 1])',...
            'linewidth',2,'color','m','marker','.')
        if exist('micronsPerPixel','var')
            plot(tempAxis(tempEventPerFrame~=inf)/framesPerSec,...
                temporalVelocity(tempEventPerFrame~=inf),'r-','linewidth',2);
            %,tempAxis(tempEventPerFrame==inf)/framesPerSec,...
            % temporalVelocity(tempEventPerFrame==inf),'bx')
            if size(hPlotNet.avVelPerEvent,2)>1
                if size(hPlotNet.avVelPerEvent,2)<6
                    title(...
                        strcat('[t_{i}, t_{f}, v] = [',...
                        (num2str([timeSlots'/framesPerSec ...
                        round(hPlotNet.avVelPerEvent')],5)),']'));
                else
                    title(...
                        strcat('[t_{i}, t_{f}, v] = [',...
                        (num2str([timeSlots(:,1:5)'/framesPerSec ...
                        round(hPlotNet.avVelPerEvent(1:5)')],5)),']'));
                end
            end
        else
            plot(tempAxis(tempEventPerFrame~=inf),...
                temporalVelocity(tempEventPerFrame~=inf),'r-o')
            title(strcat('v = [',num2str(hPlotNet.avVelPerEvent',5),']'));
        end
        
        xlabel(lab_Time_fps,'fontsize',18);
        ylabel(lab_Velo_mm_s,'fontsize',18)
    else
        plot(tempAxis(tempEventPerFrame~=inf)/framesPerSec,...
            temporalVelocity(tempEventPerFrame~=inf),'r-', ...
            'linewidth',2);
        %,tempAxis(tempEventPerFrame==inf)/framesPerSec,...
        % temporalVelocity(tempEventPerFrame==inf),'bx')
        xlabel(lab_Time_fps,'fontsize',18);
        ylabel(lab_Velo_mm_s,'fontsize',18);
        hPlotNet.avVelPerEvent=[];
    end
    axis([[0.8*min(hPlotNet.startFrame) 1.05*max(hPlotNet.stopFrame)] ...
        framesPerSec*[0.8*min(hPlotNet.avDistPerTrack) ...
        1.05*max(hPlotNet.avDistPerTrack)]]);%axis tight;
    grid on;
    
    
    
elseif (typeOfPlot==15)
    if numel(currentTrack) ~= numTracks
        selectSomeTracks = ones(numTracks,1);
    else
        selectSomeTracks = currentTrack;
    end
    
    for counterTrack=1:numTracks
        if selectSomeTracks(counterTrack)==1
            %--- determine WHICH nodes (points) correspond to the current Track
            plottingPoints = ...
                neighNetwork(1:hPlotNet.finalNetwork(counterTrack),...
                counterTrack);
            numPoints = numel(plottingPoints);
            XX = nodeNetworkX(plottingPoints);
            YY = nodeNetworkY(plottingPoints);
            ZZ = nodeNetworkZ(plottingPoints);
            
            diffBetPoints = diff([XX YY]);
            diffBetStart = [XX(1)-XX(2:end) YY(1)-YY(2:end)];
            distPerHop = sqrt((sum(diffBetPoints.^2,2)));
            distPerHop_smooth = conv(distPerHop,gaussF(13,1,1),'same');
            anglPerHop = acos(diffBetPoints(:,1)./(distPerHop+1e-30));
            distFromStart = sqrt(sum(diffBetStart.^2,2));
            
            %hPlotNet.anglFromStart(counterTrack) = atan2(...
            %(YY(end)-YY(1)) , ( (XX(end)-XX(1))+1e-30));
            hPlotNet.avDistPerTrack(counterTrack) = mean(distPerHop);
            hPlotNet.totDistPerTrack(counterTrack) = sum(distPerHop);
            %hPlotNet.maxDistPerTrack(counterTrack) = max(distPerHop);
            hPlotNet.startFrame(counterTrack) = ZZ(1);
            hPlotNet.stopFrame(counterTrack) = ZZ(end);
            %hPlotNet.avX(counterTrack) = mean(XX);
            %hPlotNet.avY(counterTrack) = mean(YY);
            %hPlotNet.anglPerHop(counterTrack) = mean(anglPerHop);
            %distV=distPerHop;
            %anglV=anglPerHop;
            %relDistP=distV(2:end)./(distV(1:end-1)+1e-30);
            %relAnglP=anglV(2:end)-anglV(1:end-1);
            %relAnglP(relDistP>1000)=[];
            %relDistP(relDistP>1000)=[];
            %hPlotNet.relDistPrevHop = [hPlotNet.relDistPrevHop;relDistP];
            %hPlotNet.relAnglPrevHop = [hPlotNet.relAnglPrevHop;relAnglP];
            
            plot3(ZZ(1:end-1),repmat(counterTrack,numPoints-1,1),...
                distPerHop_smooth,'color',colorID3(1+rem(counterTrack-1,20),:))
            hPlotNet.tempVelocity(counterTrack,round(...
                [framesPerSec*hPlotNet.startFrame(end):...
                framesPerSec*hPlotNet.stopFrame(end)]))=...
                framesPerSec*hPlotNet.avDistPerTrack(end);
            line([hPlotNet.startFrame(end)...
                hPlotNet.stopFrame(end)], ...
                framesPerSec*hPlotNet.avDistPerTrack(end)*[1 1],...
                'marker','s','linewidth',1);
        end
    end
    
    tempAxis = [framesPerSec*min(hPlotNet.startFrame):...
        framesPerSec*max(hPlotNet.stopFrame)];
    hPlotNet.tempVelocity(isnan(hPlotNet.tempVelocity))=0;
    tempVelocityPerFrame = sum(hPlotNet.tempVelocity);
    tempEventPerFrame = sum(hPlotNet.tempVelocity>0);
    tempEventPerFrame(tempEventPerFrame==0) =   inf;
    numEventsNoZero = sum(tempEventPerFrame==inf)+1;
    temporalVelocity =  tempVelocityPerFrame./tempEventPerFrame;
    q1=cumsum(temporalVelocity);
    breakPoints=find(tempEventPerFrame==inf);
    if numel(breakPoints)>0
        temporalVelocity2=temporalVelocity;
        temporalVelocity2(breakPoints)=...
            -[q1(breakPoints(1)) diff(q1(breakPoints))];
        q3=cumsum(temporalVelocity2);
        breakPoints2=[breakPoints-1 size(q3,2)];
        breakPoints3=[breakPoints2(1) -1+diff(breakPoints2)];
        
        hPlotNet.avVelPerEvent=q3(breakPoints2)./(breakPoints3+1e-100);
        timeSlots=[1 breakPoints+1;breakPoints2];
        
        timeSlots(:,(hPlotNet.avVelPerEvent==0)|...
            (hPlotNet.avVelPerEvent>1e87)|...
            (hPlotNet.avVelPerEvent<-1e87))=[];
        
        hPlotNet.avVelPerEvent(hPlotNet.avVelPerEvent==0)=[];
        hPlotNet.avVelPerEvent((hPlotNet.avVelPerEvent>1e87))=[];
        hPlotNet.avVelPerEvent((hPlotNet.avVelPerEvent<-1e87))=[];
        
        plot(timeSlots/framesPerSec,(hPlotNet.avVelPerEvent'*[1 1])',...
            'linewidth',2,'color','m','marker','.')
        if exist('micronsPerPixel','var')
            plot(tempAxis(tempEventPerFrame~=inf)/framesPerSec,...
                temporalVelocity(tempEventPerFrame~=inf),'r-','linewidth',2);
            %,tempAxis(tempEventPerFrame==inf)/framesPerSec,...
            %temporalVelocity(tempEventPerFrame==inf),'bx')
            if size(hPlotNet.avVelPerEvent,2)>1
                if size(hPlotNet.avVelPerEvent,2)<6
                    title(...
                        strcat('[t_{i}, t_{f}, v] = [',...
                        (num2str([timeSlots'/framesPerSec ...
                        round(hPlotNet.avVelPerEvent')],5)),']'));
                else
                    title(...
                        strcat('[t_{i}, t_{f}, v] = [',...
                        (num2str([timeSlots(:,1:5)'/framesPerSec ...
                        round(hPlotNet.avVelPerEvent(1:5)')],5)),']'));
                end
            end
        else
            %plot(tempAxis(tempEventPerFrame~=inf),...
            %temporalVelocity(tempEventPerFrame~=inf),'r-o')
            title(strcat('v = [',num2str(hPlotNet.avVelPerEvent',5),']'));
        end
        
        xlabel(lab_Time_fps,'fontsize',18);
        ylabel('Track ID','fontsize',18);
        zlabel(lab_Velo_mm_s,'fontsize',18)
    else
        %plot3(tempAxis(tempEventPerFrame~=inf)/framesPerSec,...
        %temporalVelocity(tempEventPerFrame~=inf),'r-','linewidth',2);
        %,tempAxis(tempEventPerFrame==inf)/framesPerSec,...
        %temporalVelocity(tempEventPerFrame==inf),'bx')
        xlabel(lab_Time_fps,'fontsize',18);
        ylabel('Track ID','fontsize',18);
        zlabel(lab_Velo_mm_s,'fontsize',18);
        hPlotNet.avVelPerEvent=[];
    end
    % axis([[0.8*min(hPlotNet.startFrame)
    % 1.05*max(hPlotNet.stopFrame)] ...
    %framesPerSec*[0.8*min(hPlotNet.avDistPerTrack) ...
    %1.05*max(hPlotNet.avDistPerTrack)]]);%axis tight;
    grid on;
    axis tight
    view(15,45);
    
    
end
