function  hPlotNet=plotTracks(handles,typeOfPlot,frameToPlot,micronsPerPixel,...
                              framesPerSec,currentColour)
%function  hPlotNet=plotTracks(handles)
%function  hPlotNet=plotTracks(handles,typeOfPlot)
%function  hPlotNet=plotTracks(handles,typeOfPlot,frameToPlot,micronsPerPixel,...
%                    framesPerSec,currentColour)
%
%--------------------------------------------------------------------------
% plotTracks  Basic command to plot the tracks of a neutrophil analysis,
%     it plots in many ways
%
%       INPUT
%         handles:          handles struct including:
%                   nodeNetwork       [numRBC detected x 12 params]
%                   finalNetwork    either 1 track or [depth of tracks x numTracks]
%
%         typeOfPlot:       a number with the following options:
%                   [-14:-1]  Will be the same options as the
%                   positive ones but with X,Y,Z
%                   [1:14]      Will plot the tracks with the
%                   following options for X,Y,T
%                                       1  = highlights longer (in distance) branches
%                                       2  = highlights faster branches
%                                       3  = highlights longer (in
%                                       number of frames) branches
%                                       4  = highlights shorter branches
%                                       5  = highlights slower branches
%                                       6  = highlights smaller branches
%                                       7  = discard branches with
%                                       small total distance,
%                                       i.e. 30% of upper half average distance
%                                       8  = discard branches with less than 3 nodes
%                                       9  = plot ONLY those
%                                       branches crossing the present Frame
%                                       10 = with labels (numbers) for the tracks
%                                       11 = all tracks in green
%                                       12 = all tracks in red
%                                       13 = top in red, bottom in green
%                                       14 = top in green, bottom in red
%                                       15 = plot [x,volume,time] as type 1
%                                       16 = plot [y,volume,time] as type 1
%                                       17 = Merge into plot
%                                       [x,volume,time] one colour per handle
%                                       18 = Merge into plot
%                                       [y,volume,time] one colour per handle
%                                       19 = Merge several separate
%                                       handles, plot one colour per handles
%
%         microsPerPixel:   the calculations of distance and velocity are made on a
%
%         framesPerPixel:   pixel/frame basis unless there is a
%         distanceProp to give proportionality
%
%       OUTPUT
%         hPlotNet: a handle to the plot with all the results
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
%     along with the PhagoSight package.  If not, see
%
%                  <http://www.gnu.org/licenses/>.
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
% method in the content,
% or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

%%

%------ no input data is received, error -------------------------
%------ at least 2 parameters are required
if (nargin <1)
    help plotTracks;
    hPlotNet=[];
    return;
end
%------ arguments received revision   ----------------------------
if (~exist('typeOfPlot','var'))
    typeOfPlot=1;
end
if (~exist('frameToPlot','var') || isempty(frameToPlot))
    frameToPlot=1;
end
if (~exist('currentColour','var'));
    currentColour=1;
end
changeSecToMin = 0;

if isa(handles,'cell')
    cla;
    [numHRows,numHCols] = size(handles);
    for counterHandles =1:max(numHRows,numHCols)
        if frameToPlot==1
            plotTracks(handles{counterHandles},typeOfPlot,...
                       frameToPlot,1,1,counterHandles);
        else
            plotTracks(handles{counterHandles},typeOfPlot,...
                       handles{counterHandles}.distanceNetwork.numHops>...
                       frameToPlot,1,1,counterHandles);
        end
    end
    axis tight;
else

    typeOfPlot3 = 0;
    if (typeOfPlot<0)
        nodeNetwork = handles.nodeNetwork(:,[1 2 3]);
        typeOfPlot = abs(typeOfPlot);
        typeOfPlot2 = 1;
    else
        if (typeOfPlot >= 15 && typeOfPlot <= 18) %Replaced with equivalent condition
            switch typeOfPlot
              case 15
                nodeNetwork = handles.nodeNetwork(:,[2 10 5 ]);
                typeOfPlot3 = 1;
                typeOfPlot = 1;
                typeOfPlot2 = 0;

              case 16
                nodeNetwork = handles.nodeNetwork(:,[1 10 5 ]);
                typeOfPlot3 = 2;
                typeOfPlot = 1;
                typeOfPlot2 = 0;

              case 17
                nodeNetwork = handles.nodeNetwork(:,[2 10 5 ]);
                typeOfPlot3 = 1;
                clear micronsPerPixel framesPerSec;
              case 18
                nodeNetwork = handles.nodeNetwork(:,[1 10 5 ]);
                typeOfPlot3 = 2;
                clear micronsPerPixel framesPerSec;

            end
            typeOfPlot2 = 0;

        else
            nodeNetwork = handles.nodeNetwork(:,[1 2 5]);
            typeOfPlot2 = 0;
        end
    end

    % the number of tracks can be restricted to a subset of them
    % with certain characteristics by addressing the matrix for
    % instance:
    % plotTracks(handles,15,(handles.distanceNetwork.numHops>30))
    % thus frameToPlot acts as the definition of the subset EXCEPT if the type of
    % plot is 9 and in that case they act as the coordinates in which the tracks
    % start and stop
    if typeOfPlot==9
        finalNetwork = handles.finalNetwork;
    else
        if numel(frameToPlot)>1

            if any(frameToPlot)
                try
                    finalNetwork = handles.finalNetwork(:,frameToPlot);
                catch
                    finalNetwork = handles.finalNetwork;
                end
                frameToPlot = 1;
            else
                %the conditions have removed all the tracks from the plot, end
                hPlotNet = [];
                return;
            end
        else
            finalNetwork = handles.finalNetwork;
        end
    end

    %define all the labels according to the input arguments
    if (~exist('micronsPerPixel','var'))||(~exist('framesPerSec','var'))
        lab_Velo_mm_s  = 'Velocity [pix/frame]';
    else
        if framesPerSec<0.05

            lab_Velo_mm_s = 'Velocity [\mum /min]';
        else
            lab_Velo_mm_s = 'Velocity [\mum /s]';
        end
    end

    if ~exist('micronsPerPixel','var') ;
        micronsPerPixel = 1;
        if typeOfPlot3>0
            lab_Cols_dist = 'Volume [voxels]';
            switch typeOfPlot3
              case 1
                nodeNetwork = handles.nodeNetwork(:,[2 10 5 ]);
                lab_Rows_dist = 'Columns [pixels]';
              case 2
                nodeNetwork = handles.nodeNetwork(:,[1 10 5 ]);
                lab_Rows_dist = 'Rows [pixels]';
            end
        else
            lab_Cols_dist = 'Columns [pixels]';
            lab_Rows_dist = 'Rows [pixels]';
        end
        yLabUnits =' [pixels]';
    else
        lab_Rows_dist = 'Rows [ \mum ]';
        lab_Cols_dist = 'Columns [ \mum ]';
        yLabUnits =' [\mu m]';
    end
    if (typeOfPlot2==0)
        if ~exist('framesPerSec','var') ;
            framesPerSec = 1;
            lab_Time_fps = 'Time [frames]'  ;
        elseif framesPerSec<0.05
            changeSecToMin = 1;
            lab_Time_fps = 'Time [min]';
        else
            lab_Time_fps = 'Time [sec]';
        end
    else
        lab_Time_fps = 'Z-Position [slices]' ;
        if ~exist('framesPerSec','var') ;
            framesPerSec = 1;
        end
    end

    if changeSecToMin==1
        framesPerSec = framesPerSec * 60;
    end

    %------ regular size check ---------------------------------------
    numTracks = size(finalNetwork,2);

    % Get axis ready
    h0 = gcf;
    h1 = gca;
    %cla;
    hChild = allchild(h0);
    hColP = findobj(hChild,'tag','colorBarPlot');
    hColB = findobj(hChild,'tag','colorbarLabel');
    set(hColP,'visible','off');
    set(hColB,'visible','off');
    hold on;
    set(h1,'Color',[1 1 1]);


    %------- colour codes for the plot
    colorID2 = {'o', 'x', 's', '+', 'v',  '*', '^', '<', '>','h','p'};
    colorID3=0.9*[ [0 0 1];[0 1 0];[1 0 0];...
                   [1 0 0.75]; [0 1 0.25]; [0 0 0];...
                   [0.2 1 0.2];    [0.5 0.75 0.5]; [1 0.5 1];...
                   [1 0 0];    [0.75 0 0]; [0.35 0 0];...
                   [0 1 1];    [1 0 1];    [1 1 0 ];...
                   [0.5 0.5 0];[0.5 0 0.5];[0 0.5 0.5];...
                   [0.25 0.25 0.25];[0.25 1 0];[0.25 0 0.25];[0 0.25 0.25]];


    %-------- The structure of nodeNetwork will be [X Y z] each
    %-------- column is separately transformed int
    %-------- a matrix [depthB x numTracks] which will be used as
    %X,Y,Z for the plot3 function
    %-------- how to plot some nets that DO NOT have same number of nodes????
    nodeNetworkX = (nodeNetwork(:,2));
    nodeNetworkY = (nodeNetwork(:,1));
    nodeNetworkZ = (nodeNetwork(:,3));


    % ------------------ Calibration of the measurements -------------------------
    nodeNetworkZ = nodeNetworkZ/framesPerSec;
    nodeNetworkX = nodeNetworkX * micronsPerPixel;
    nodeNetworkY = nodeNetworkY * micronsPerPixel;

    %----- a neighbourhood for the group of tracks has been received
    neighNetwork = finalNetwork;
    %----- finalNetwork will keep a vector with the lengths of the branches
    hPlotNet.finalNetwork = sum(neighNetwork>0);
    hPlotNet.numTracks = numTracks;
    hPlotNet.relDistPrevHop = [];
    hPlotNet.relAnglPrevHop = [];

    %--------- loop over the tracks -----------------------------------
    for counterTrack=1:numTracks
        %------ determine WHICH nodes (points) correspond to the current Track
        plottingPoints = ...
            neighNetwork(1:hPlotNet.finalNetwork(counterTrack),counterTrack);
        if numel(plottingPoints)>1
            if typeOfPlot3>0
                tempVec = repmat(nodeNetworkX(plottingPoints),[1,3]);
                tempVec(1:end-1,2) = nodeNetworkX(plottingPoints(2:end));
                tempVec(2:end,3) = nodeNetworkX(plottingPoints(1:end-1));
                XX = mean(tempVec,2);
            else
                XX = nodeNetworkX(plottingPoints);
            end

            YY = nodeNetworkY(plottingPoints);
            ZZ = nodeNetworkZ(plottingPoints);
            ZZ2 = handles.nodeNetwork(plottingPoints,3);
            TT = handles.nodeNetwork(plottingPoints,5);
            diffBetPoints = diff([XX YY]);
            diffBetStart = [XX(1)-XX(2:end) YY(1)-YY(2:end)];
            distPerHop = sqrt((sum(diffBetPoints.^2,2)));
            anglPerHop = acos(diffBetPoints(:,1)./(distPerHop+1e-30));
            distFromStart = sqrt(sum(diffBetStart.^2,2));

            hPlotNet.anglFromStart(counterTrack) = atan2( (YY(end)-YY(1)) , ...
                                                          ( (XX(end)-XX(1))+1e-30));
            hPlotNet.avDistPerTrack(counterTrack) = mean(distPerHop);
            hPlotNet.totDistPerTrack(counterTrack) = sum(distPerHop);
            hPlotNet.maxDistPerTrack(counterTrack) = max(distPerHop);
            hPlotNet.startFrame(counterTrack) = ZZ(1);
            hPlotNet.stopFrame(counterTrack) = ZZ(end);
            hPlotNet.startX(counterTrack) = XX(1);
            hPlotNet.stopX(counterTrack) = XX(end);
            hPlotNet.startY(counterTrack) = YY(1);
            hPlotNet.stopY(counterTrack) = YY(end);

            hPlotNet.startFrameT(counterTrack) = TT(1);
            hPlotNet.stopFrameT(counterTrack) = TT(end);

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
            hPlotNet.relAnglPrevHop =[hPlotNet.relAnglPrevHop;relAnglP];

            if typeOfPlot<11
                hPlotNet.handlePlot(counterTrack) = plot3(...
                    XX,YY,ZZ,'marker',colorID2{1+rem(counterTrack-1,11)},...
                    'color',colorID3(1+rem(counterTrack-1,20),:),'markersize',...
                    4,'linewidth',1);
                if typeOfPlot==10
                    %------ add a label with the number of the track
                    text(2*XX(1)-XX(2)+ceil(handles.rows/100),...
                         2*YY(1)-YY(2)+ceil(handles.cols/100),...
                         0.9*ZZ(1),num2str(counterTrack));
                end
            elseif typeOfPlot==11
                % all tracks green
                hPlotNet.handlePlot(counterTrack) = plot3(...
                    XX,YY,ZZ,'marker','.','color','g','markersize',4,'linewidth',2);
            elseif typeOfPlot==12
                % all tracks red
                hPlotNet.handlePlot(counterTrack) = plot3(...
                    XX,YY,ZZ,'marker','.','color','r','markersize',4,'linewidth',2);
            elseif typeOfPlot==13
                %Distinguish between top and bottom according to
                %the ChannelDistribution
                if ~isfield(handles,'ChannelDistribution')
                    hPlotNet.handlePlot(counterTrack) = plot3(...
                        XX,YY,ZZ,'marker','.','color','g','markersize',...
                        4,'linewidth',2);
                else

                    if (mean(ZZ2)<=(handles.ChannelDistribution(2)))&&...
                            (mean(ZZ2)>=(handles.ChannelDistribution(1)))
                        hPlotNet.handlePlot(counterTrack) = plot3(...
                            XX,YY,ZZ,'marker','.','color','g','markersize',...
                            4,'linewidth',2);
                    else
                        hPlotNet.handlePlot(counterTrack) = plot3(...
                            XX,YY,ZZ,'marker','.','color','r','markersize',...
                            4,'linewidth',2);
                    end
                end
            elseif typeOfPlot==14
                if ~isfield(handles,'ChannelDistribution')
                    hPlotNet.handlePlot(counterTrack) = plot3(...
                        XX,YY,ZZ,'marker','.','color','g','markersize',...
                        4,'linewidth',2);
                else

                    if (mean(ZZ2)<=(handles.ChannelDistribution(2)))...
                            &&(mean(ZZ2)>=(handles.ChannelDistribution(1)))
                        hPlotNet.handlePlot(counterTrack) = plot3(...
                            XX,YY,ZZ,'marker','.','color','r','markersize',...
                            4,'linewidth',2);
                    else
                        hPlotNet.handlePlot(counterTrack) = plot3(...
                            XX,YY,ZZ,'marker','.','color','g','markersize',...
                            4,'linewidth',2);
                    end
                end
            elseif typeOfPlot>=17
                % all tracks red
                hPlotNet.handlePlot(counterTrack) = plot3(...
                    XX,YY,ZZ,'marker','.','color',...
                    colorID3(1+rem(currentColour-1,20),:),...
                    'markersize',4,'linewidth',2);
            end
        end
    end


    zlabel(lab_Time_fps,'fontsize',18);
    xlabel(lab_Cols_dist,'fontsize',15);
    ylabel(lab_Rows_dist,'fontsize',15);

    hPlotNet.totDistPerTrack(isnan(hPlotNet.totDistPerTrack))=0;
    hPlotNet.avDistPerTrack(isnan(hPlotNet.avDistPerTrack))=0;
    hPlotNet.finalNetwork(isnan(hPlotNet.finalNetwork))=0;

    if typeOfPlot<10
        colorID6 = interp1(jet,linspace(1,64,hPlotNet.numTracks));

        widthID = (linspace(2.5,2.5,hPlotNet.numTracks));
        switch typeOfPlot
          case 1
            %------- highlights longer (in distance) branches
            [m1,m2] = sort(hPlotNet.totDistPerTrack);
            initBranch = 1;
          case 2
            %------- highlights faster branches    -----------------------------
            [m1,m2] = sort(hPlotNet.avDistPerTrack);
            initBranch = 1;
            if exist('framesPerSec','var')
                colorbar('tag','colorBarPlot','ytick',(1:9:64),...
                          'yticklabel',num2str(...
                              framesPerSec*(linspace(m1(1),m1(end),8))',4));
                annotation1 = annotation(h0,'textbox','tag','colorbarLabel',...
                                         'Position',[0.83 0.03 0.09 0.07],...
                                         'String',{lab_Velo_mm_s},'linestyle',...
                                         'none','fontsize',12,...
                                         'FitHeightToText','on');
            else
                colorbar('tag','colorBarPlot','ytick',(1:9:64),...
                         'yticklabel',num2str((linspace(m1(1),m1(end),8))',3));
                annotation1 = annotation(h0,'textbox','tag','colorbarLabel',...
                                         'Position',[0.845 0.02 0.095 0.07],...
                                         'String',{lab_Velo_mm_s},'linestyle',...
                                         'none','fontsize',14,...
                                         'FitHeightToText','on');

            end
            colormap jet;

          case 3
            %------- highlights longer (in number of frames) branches
            [m1,m2] = sort(hPlotNet.finalNetwork);
            initBranch = 1;
                case 4
                  %------- highlights shorter branches
                  [m1,m2] = sort(hPlotNet.totDistPerTrack);
                  m2 = m2(end:-1:1);
                  initBranch = 1;
          case 5
            %------- highlights slower branches
            [m1,m2] = sort(hPlotNet.avDistPerTrack);
            m2 = m2(end:-1:1);
            initBranch = 1;
          case 6
            %------- highlights smaller branches
            [m1,m2] = sort(hPlotNet.finalNetwork);
            m2 = m2(end:-1:1);
            initBranch = 1;
          case 7
            %------- discard branches with small total distance,
            %i.e. 30% of upper half average distance
            [m1,m2] = sort(hPlotNet.totDistPerTrack);
            initBranch = max(1,sum(m1<floor((ceil(mean(m1(ceil(end/2):end))))*0.3)));
            widthID = (0.5+logspace(0.5,4.8,hPlotNet.numTracks)/30100);
          case 8
            %------- discard branches with less than 3 nodes
            [m1,m2] = sort(hPlotNet.finalNetwork);
            initBranch = max(1,sum(m1<4));
            widthID = (0.2+logspace(0.5,4.8,hPlotNet.numTracks)/30100);
          case 9

            if numel(frameToPlot)==4
                tracksToKeepX = find(...
                    (hPlotNet.startX>=frameToPlot(1))&...
                    (hPlotNet.stopX<=frameToPlot(2)));
                tracksToKeepY = find(...
                    (hPlotNet.startY>=frameToPlot(3))&...
                    (hPlotNet.stopY<=frameToPlot(4)));

                tracksToKeep = intersect(tracksToKeepX,tracksToKeepY);
                tracksNotToKeep = setdiff((1:numTracks),tracksToKeep);
                m2 = [tracksToKeep tracksNotToKeep];

                widthID(numel(tracksToKeep)+1:numTracks)=0.5;
                colorID6(numel(tracksToKeep)+1:numTracks,:)=...
                    repmat([1 1 1]*0.95,[numTracks-numel(tracksToKeep),1]);

            else

                %------- plot ONLY those branches crossing the present Frame
                tracksToKeep = find(...
                    (hPlotNet.startFrameT<=frameToPlot(1))&...
                    (hPlotNet.stopFrameT>=frameToPlot(1)));
                tracksToKeep2 = find(...
                    (hPlotNet.startFrameT<=(frameToPlot(1)+1))&...
                    (hPlotNet.stopFrameT>=(frameToPlot(1)+1)));
                tracksToKeep3 = setdiff(tracksToKeep2,tracksToKeep);
                tracksNotToKeep1 = setdiff((1:numTracks),tracksToKeep);
                tracksNotToKeep2 = setdiff(tracksNotToKeep1,tracksToKeep3);

                m2 = [tracksToKeep3 tracksNotToKeep2 tracksToKeep];

                widthID(numel(tracksToKeep3)+1:numTracks- ...
                        numel(tracksToKeep)) = 0.5;
                colorID6(numel(tracksToKeep3)+1:numTracks- ...
                         numel(tracksToKeep),:) = ...
                    repmat([1 1 1]*0.8,[numTracks-...
                                    numel(tracksToKeep)-...
                                    numel(tracksToKeep3),1]);

            end

            initBranch =1;
        end
        hPlotNet.sortedBranches = m2;

        if hPlotNet.numTracks > 15
            sizeForMarker = 5;
        else
            sizeForMarker = 9;
        end

        for counterTrack=initBranch:hPlotNet.numTracks
            if (hPlotNet.handlePlot(m2(counterTrack)))~=0
                set(hPlotNet.handlePlot(m2(counterTrack)),...
                    'linewidth',widthID(counterTrack),'color',...
                    colorID6(counterTrack,:),'marker','.');
                if typeOfPlot==1
                    plot3(hPlotNet.stopX(m2(counterTrack)),...
                          hPlotNet.stopY(m2(counterTrack)),...
                          hPlotNet.stopFrame(m2(counterTrack)),...
                          '^','color','r','markerfacecolor',...
                          colorID6(counterTrack,:),'markersize',sizeForMarker)
                    plot3(hPlotNet.startX(m2(counterTrack)),...
                          hPlotNet.startY(m2(counterTrack)),...
                          hPlotNet.startFrame(m2(counterTrack)),...
                          's','color','b','markerfacecolor',...
                          colorID6(counterTrack,:),'markersize',sizeForMarker-1)
                end
            end
        end
    end

    grid on;rotate3d on;
    axis ij
    view(20,50)
    if typeOfPlot3>0
        axis([1 max(handles.nodeNetwork(:,10)) 1 ...
              handles.rows*micronsPerPixel 1 handles.numFrames/framesPerSec])
    else
        if typeOfPlot2==0
            axis([1 handles.cols*micronsPerPixel 1 ...
                  handles.rows*micronsPerPixel 1 handles.numFrames/framesPerSec])
        else
            axis([1 handles.cols*micronsPerPixel 1 ...
                  handles.rows*micronsPerPixel 1 handles.levs])
        end
    end
end
