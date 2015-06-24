function [handles] = verifyThresholdNeutrophils(handles,dataR,indexR,indexC,indexL,askUser)
%function [handles] = verifyThresholdNeutrophils(handles,dataR,indexR,indexC,indexL,askUser)
%
%--------------------------------------------------------------------------
% verifyThresholdNeutrophils  is used to display the datasets with the
%     currently selected thresholds for the user to verify those ones.
%
%       INPUT
%           handles:	contains the current thresholds and dimensions of data
%           dataR:      the data to be thresholded 
%           indexR:     the region of interest (in Rows) to display a zoom
%           indexC:     the region of interest (in Columns) to display a zoom
%           indexL:     the region of interest (in Levels) to display a zoom
%                       these are calculated to show the region of max intensity
%
%       OUTPUT
%           handles:    updated handles struct with the thresholds re-calculated 
%                       and the handles to the figure (to close them later)
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
%find the position where the maximum intensity occurs

indexLevels = [handles.ChannelDistribution(1):handles.ChannelDistribution(2) ...
    handles.ChannelDistribution(5):handles.ChannelDistribution(6)];

indexLevels (indexLevels==0) = [];

indexGreen = handles.ChannelDistribution(1):handles.ChannelDistribution(2);
indexRed = handles.ChannelDistribution(5):handles.ChannelDistribution(6);
[rows,cols,levs] = size(dataR);
dataIn = dataR;

if numel(handles.thresLevel)==2
%%
    maxIntensityProj                                        = max(dataR(:,:,indexLevels),[],3);
    dataR_LowThres                                          = maxIntensityProj>handles.thresLevel(1);
    dataR_HighThres                                         = maxIntensityProj>handles.thresLevel(2);
    % calculate the size of the neutrophils so far to fix min blob
    lowT                                                    = min(handles.thresLevel);
    highT                                                   = max(handles.thresLevel);
    GreenLowT                                               = bwlabeln(dataR(:,:,indexLevels)>lowT);
    neutropToKeep                                           = unique(GreenLowT.*(dataR(:,:,indexLevels)>highT));
    [dataL_green,numNeutropGreen]                           = bwlabeln(ismember(GreenLowT,neutropToKeep(2:end))); %#ok<NASGU>
    
    regS                                                    = regionprops(dataL_green);
    neutsVolumes                                            = sort(round([regS.Area]));
    neutsVolumes                                            = neutsVolumes (end:-1:1);
    maxNum2show                                             = min(numNeutropGreen,8);
    
    dataR_2Levels                                           = dataR_LowThres+dataR_HighThres;
    
%%    
else
    %the only way to have 4 threshold levels is: 
    %a) once they have been determined here, or b) they were provided as input 
    maxIntensityProj_G                                      = max(dataR(:,:,indexGreen),[],3);
    dataR_LowThres_G                                        = maxIntensityProj_G>handles.thresLevel(1);
    dataR_HighThres_G                                       = maxIntensityProj_G>handles.thresLevel(2);
    
    maxIntensityProj_R                                      = max(dataR(:,:,indexRed),[],3);
    dataR_LowThres_R                                        = maxIntensityProj_R>handles.thresLevel(3);
    dataR_HighThres_R                                       = maxIntensityProj_R>handles.thresLevel(4);

    
        % calculate the size of the neutrophils so far to fix min blob
    lowT                                                    = min(handles.thresLevel(1:2));
    highT                                                   = max(handles.thresLevel(1:2));
    GreenLowT                                               = bwlabeln(dataR(:,:,indexGreen)>lowT);
    neutropToKeep                                           = unique(GreenLowT.*(dataR(:,:,indexGreen)>highT));
    [dataL_green,numNeutropGreen]                           = bwlabeln(ismember(GreenLowT,neutropToKeep(2:end))); 
    
    regS                                                    = regionprops(dataL_green);
    neutsVolumes_G                                          = sort(round([regS.Area]));
    neutsVolumes_G                                          = neutsVolumes_G (end:-1:1);
   % maxNum2show                                             = min(numNeutropGreen,8);

    lowT                                                    = min(handles.thresLevel(3:4));
    highT                                                   = max(handles.thresLevel(3:4));
    RedLowT                                                 = bwlabeln(dataR(:,:,indexRed)>lowT);
    neutropToKeepR                                          = unique(RedLowT.*(dataR(:,:,indexRed)>highT));
    [dataL_Red,numNeutropRed]                               = bwlabeln(ismember(RedLowT,neutropToKeepR(2:end))); 

    regS_R                                                  = regionprops(dataL_Red);
    neutsVolumes_R                                          = sort(round([regS_R.Area]));
    neutsVolumes_R                                          = neutsVolumes_R (end:-1:1);
    neutsVolumes                                            = [neutsVolumes_R neutsVolumes_G];
    
    maxNum2show                                             = min(numNeutropRed+numNeutropGreen,8);
    dataR_2Levels_temp                                           = dataR_LowThres_G+dataR_HighThres_G +(3*dataR_LowThres_R+dataR_HighThres_R);
    dataR_2Levels(:,:,1) = (dataR_LowThres_R-dataR_HighThres_R) + (dataR_LowThres_G-dataR_HighThres_G) + dataR_HighThres_R;
    dataR_2Levels(:,:,2) = (dataR_LowThres_R-dataR_HighThres_R) + (dataR_LowThres_G-dataR_HighThres_G) + dataR_HighThres_G;
    dataR_2Levels(:,:,3) = (dataR_LowThres_R-dataR_HighThres_R) + (dataR_LowThres_G-dataR_HighThres_G);
    
    dataR_2Levels(dataR_2Levels>1) = 1;
    
    maxIntensityProj                                        = maxIntensityProj_G +maxIntensityProj_R;

    
end
dataR                                                       = dataR(:,:, indexLevels);

%%
%define a region to display
maxDataR = floor(max(dataR(:)));
minDataR = floor(min(dataR(:)));
maxIndex = find(dataR>=maxDataR);
[rr,cc,ll] = ind2sub([rows, cols, levs],maxIndex);

if (exist('indexR','var') && exist('indexC','var') && exist('indexL','var') &&...
        numel(indexR)==0 && numel(indexC)==0 && numel(indexL)==0)
   %These means the arguments are empty and they were used to provide a
   %value for the argument askUser. All three variable must be delete from
   %the workspace to allow this function running normally.
   clear indexR;
   clear indexC;
   clear indexL
end

if ~exist('indexR','var')
    indexR                                                  = max(1,rr-50):min(rows,rr+50);
    indexC                                                  = max(1,cc-50):min(cols,cc+50);
    indexL                                                  = 1:size(dataR,3);
end
%%
if isempty(indexR);  indexR = 1 : 10; end
if isempty(indexC);  indexC = 1 : 10; end

%%
fig1 = figure(31);
clf
jet4=jet;
set(fig1,'userdata',jet4)
set(fig1,'renderer','opengl')
set(fig1,'position',[60,280,1300,600])
set(fig1,'menuBar','none','toolBar','none','numbertitle','off')

%% Compute the histograms for three points in time, first, halfway and last and display

dir1 = dir(strcat(handles.dataRe,'/*.mat'));
timePointsToRead = round(linspace(1,size(dir1,1),3));

dataRMid = load(strcat(handles.dataRe,'/',...
        dir1(timePointsToRead(2)).name),'dataR');
dataRMid = double(dataRMid.dataR);

dataRFin = load (strcat(handles.dataRe,'/',...
        dir1(timePointsToRead(3)).name),'dataR');
dataRFin = double(dataRFin.dataR);

%% get the histograms according to green / red
if (indexGreen(1)==0)||(indexRed(1)==0)
    %only one channel present, process as a single case
    dataRMid = dataRMid(:,:,indexLevels);
    dataRFin = dataRFin(:,:,indexLevels);
    minDataRM = floor(min(dataRMid(:)));
    maxDataRM = floor(max(dataRMid(:)));
    minDataRF = floor(min(dataRFin(:)));
    maxDataRF = floor(max(dataRFin(:)));
    
    [yHistI,xHistI] = hist(dataR(:)   ,linspace(minDataR,maxDataR,100));
    [yHistM,xHistM] = hist(dataRMid(:),linspace(minDataRM,maxDataRM,100));
    [yHistF,xHistF] = hist(dataRFin(:),linspace(minDataRF,maxDataRF,100));
else
    %both green and red are present, process separately
    
    dataRMid_G = dataRMid(:,:,indexGreen);
    dataRFin_G = dataRFin(:,:,indexGreen);
    dataRIni_G = dataIn(:,:,indexGreen);
    maxDataRM_G = floor(max(dataRMid_G(:)));
    maxDataRF_G = floor(max(dataRFin_G(:)));
    minDataRM_G = floor(min(dataRMid_G(:)));
    minDataRF_G = floor(min(dataRFin_G(:)));
  
    [yHistI_G,xHistI_G] = hist(dataRIni_G(:),linspace(minDataR,maxDataR,100));
    [yHistM_G,xHistM_G] = hist(dataRMid_G(:),linspace(minDataRM_G,maxDataRM_G,100));
    [yHistF_G,xHistF_G] = hist(dataRFin_G(:),linspace(minDataRF_G,maxDataRF_G,100));

    dataRMid_R = dataRMid(:,:,indexRed);
    dataRFin_R = dataRFin(:,:,indexRed);
    dataRIni_R = dataIn(:,:,indexRed);
    maxDataRM_R = floor(max(dataRMid_R(:)));
    maxDataRF_R = floor(max(dataRFin_R(:)));
    minDataRM_R = floor(min(dataRMid_R(:)));
    minDataRF_R = floor(min(dataRFin_R(:)));
    
    [yHistI_R,xHistI_R] = hist(dataRIni_R(:),linspace(minDataR,maxDataR,100));
    [yHistM_R,xHistM_R] = hist(dataRMid_R(:),linspace(minDataRM_R,maxDataRM_R,100));
    [yHistF_R,xHistF_R] = hist(dataRFin_R(:),linspace(minDataRF_R,maxDataRF_R,100));
end

%% Draw the labelled and intensity images plus ROI

h231=subplot(231);set(h231,'position',[ 0.03 0.535     0.270    0.4000]);
imagesc(dataR_2Levels)
drawsquare_r(indexC(1), indexR(1),indexC(end), indexR(end))
title(strcat('(a) Volume Neuts (top ',num2str(maxNum2show),'/',...
    num2str(numNeutropGreen),'):',num2str(neutsVolumes(1:maxNum2show))),...
    'fontsize',10)

h232=subplot(232);set(h232,'position',[ 0.33 0.535     0.270    0.4000]);
imagesc(maxIntensityProj)
drawsquare_r(indexC(1), indexR(1),indexC(end), indexR(end))
colorbar;
title('(b)','fontsize',18)

h234=subplot(234);set(h234,'position',[ 0.03 0.045    0.270    0.4000]);
imagesc(dataR_2Levels(indexR,indexC,:))
title('(c)','fontsize',18)

h235=subplot(235);set(h235,'position',[ 0.33 0.045    0.270    0.4000]);
imagesc(maxIntensityProj(indexR,indexC))
colorbar;
title('(d)','fontsize',18)



h233=subplot(233);set(h233,'position',[ 0.65 0.6    0.310    0.32]);


%% arrange the plots of the histograms according to green / red channels

if (indexGreen(1)==0)||(indexRed(1)==0)

    xAxisMin = min([xHistI xHistM xHistF]);
    xAxisMax = max([xHistI xHistM xHistF]);

    hold off
    plot(xHistI,log10((yHistI+1)/sum(yHistI)),'b-','linewidth',2)
    hold on
    plot(xHistM,log10((yHistM+1)/sum(yHistM)),'m--','linewidth',1)
    plot(xHistF,log10((yHistF+1)/sum(yHistF)),'k-.','linewidth',1)
    legend('Time frame = 1',strcat('Time frame = ',num2str(timePointsToRead(2))),strcat('Time frame = ',num2str(timePointsToRead(3))))
    
    
    hold on
    plot([handles.thresLevel(1) handles.thresLevel(1)],[-18 0],'r')
    plot([handles.thresLevel(2) handles.thresLevel(2)],[-18 0],'r')
    axis tight;
    grid on
    ylabel('log10 of the occurrence')
    axis ([xAxisMin xAxisMax -0.5+log10(min((1+yHistI)/sum(yHistI))) 0])
    title('(e)','fontsize',18)
else
    
    xAxisMin = min([xHistI_G xHistM_G xHistF_G xHistI_R xHistM_R xHistF_R]);
    xAxisMax = max([xHistI_G xHistM_G xHistF_G xHistI_R xHistM_R xHistF_R]);

    hold off
    plot(xHistI_G,log10((yHistI_G+1)/sum(yHistI_G)),...
            '-' ,'linewidth',2,'color',[0 0.61 0])
    hold on
    plot(xHistM_G,log10((yHistM_G+1)/sum(yHistM_G)),...
        '--','linewidth',1,'color',[0 0.51 0])
    plot(xHistF_G,log10((yHistF_G+1)/sum(yHistF_G)),...
        '-.','linewidth',1,'color',[0 0.41 0])

    plot(xHistI_R,log10((yHistI_R+1)/sum(yHistI_R)),...
        '-' ,'linewidth',2,'color',[0.91 0 0])
    plot(xHistM_R,log10((yHistM_R+1)/sum(yHistM_R)),...
        '--','linewidth',1,'color',[0.81 0 0])
    plot(xHistF_R,log10((yHistF_R+1)/sum(yHistF_R)),...
        '-.','linewidth',1,'color',[0.71 0 0])

    
    hLegend = legend('Green: Time frame = 1',...
            strcat('Green: Time frame = ',num2str(timePointsToRead(2))),...
            strcat('Green: Time frame = ',num2str(timePointsToRead(3))),...
            'Red: Time frame = 1',...
            strcat('Red: Time frame = ',num2str(timePointsToRead(2))),...
            strcat('Red: Time frame = ',num2str(timePointsToRead(3))));
    set(hLegend,'fontsize',9);
    
    if numel(handles.thresLevel)==2
        plot([handles.thresLevel(1) handles.thresLevel(1)],[-18 0],'b')
        plot([handles.thresLevel(2) handles.thresLevel(2)],[-18 0],'b')
    else
        plot([handles.thresLevel(1) handles.thresLevel(1)],[-18 0],'g')
        plot([handles.thresLevel(2) handles.thresLevel(2)],[-18 0],'g')
        plot([handles.thresLevel(3) handles.thresLevel(3)],[-18 0],'r')
        plot([handles.thresLevel(4) handles.thresLevel(4)],[-18 0],'r')
        
    end
    axis tight;
    grid on
    ylabel('log10 of the occurrence')
    axis ([xAxisMin xAxisMax -0.5+log10(min((1+yHistI_G)/sum(yHistI_G))) 0])
    title('(e)','fontsize',18)
    
    
end

%%
if numel(indexL)>1
    h236=subplot(236);set(h236,'position',[ 0.64 0.045    0.320    0.45000]);
    cla
    p1 = patch(isosurface(dataR(indexR,indexC,indexL), handles.thresLevel(1)),'FaceColor','c','EdgeColor','none');
    p2 = patch(isosurface(dataR(indexR,indexC,indexL), handles.thresLevel(2)),'FaceColor','m','EdgeColor','none');
    
    axis ij
    view(10,70)
    alpha (0.5)
    q1=camlight('right');
    set(q1,'Position',[-320,-260,0])
    q2=camlight('right');
    set(q2,'Position',[540,-100,60])
    %camlight left; camlight; 
    lighting gouraud
    rotate3d on
    grid on
    title('(f)','fontsize',18)

end


%% set the colour maps
h_uimenu2 = uimenu(gcf,'Label','Colour Map');
h_uimenu2_jet       =   uimenu(h_uimenu2,'Label','jet',       'Callback','colormap(jet);');

h_uimenu2_jet2      =   uimenu(h_uimenu2,'Label','jet (low)', 'Callback','qq(:,3) =[linspace(0.5,1,8)'';ones(16,1);linspace(1,0,15)'';zeros(25,1)] ;qq(:,1) = [zeros(24,1);linspace(0,1,15)'';ones(16,1);linspace(1,0.5,9)''];qq(:,2) = [zeros(8,1);linspace(0,1,15)'';ones(16,1);linspace(1,0,15)'';zeros(10,1)]; colormap(qq([1:3:40 41:end],:));clear qq;');
h_uimenu2_jet3      =   uimenu(h_uimenu2,'Label','jet (high)','Callback','qq(:,3) =[linspace(0.5,1,8)'';ones(16,1);linspace(1,0,15)'';zeros(25,1)] ;qq(:,1) = [zeros(24,1);linspace(0,1,15)'';ones(16,1);linspace(1,0.5,9)''];qq(:,2) = [zeros(8,1);linspace(0,1,15)'';ones(16,1);linspace(1,0,15)'';zeros(10,1)]; colormap(qq([1:20 21:2:end],:));clear qq;');
h_uimenu2_hot       =   uimenu(h_uimenu2,'Label','hot',       'Callback','colormap(hot);');
h_uimenu2_hot2      =   uimenu(h_uimenu2,'Label','hot^2',     'Callback','colormap(hot.^2);');
h_uimenu2_hot12     =   uimenu(h_uimenu2,'Label','hot^(1/2)', 'Callback','colormap(hot.^(0.5));');
h_uimenu2_gray      =   uimenu(h_uimenu2,'Label','gray',      'Callback','colormap(gray);');
h_uimenu2_gray2     =   uimenu(h_uimenu2,'Label','gray^2',    'Callback','colormap(gray.^2);');
h_uimenu2_gray12    =   uimenu(h_uimenu2,'Label','gray^(1/2)','Callback','colormap(gray.^(0.5));');
h_uimenu2_green     =   uimenu(h_uimenu2,'Label','green','Callback','qq=hot; qq2(:,2)=qq(:,1); qq2(:,3)=qq(:,3);  qq2(:,1)=qq(:,3); colormap(qq2);');
h_uimenu2_green2    =   uimenu(h_uimenu2,'Label','green^2','Callback','qq=hot; qq2(:,2)=qq(:,1); qq2(:,3)=qq(:,3);  qq2(:,1)=qq(:,3); colormap(qq2.^2);');
h_uimenu2_green12   =   uimenu(h_uimenu2,'Label','green^(1/2)','Callback','qq=hot; qq2(:,2)=qq(:,1); qq2(:,3)=qq(:,3);  qq2(:,1)=qq(:,3); colormap(qq2.^(0.5));');
                    

%These conditions were added to avoid confirmation from the user, this way is
%easier to test the method
if(~exist('askUser','var'))
   askUser = true; 
end
if (askUser)

    reply = input('Do you want modify plot or threshold levels? Y/N [n]: ','s');
    if (strcmp(reply,'y'))||(strcmp(reply,'Y'))
        disp ('Type the new values for the following parameter or press RETURN to leave unchanged');

        %determine if separate thresholds will be used for green/red
        if (indexGreen(1)~=0)&&(indexRed(1)~=0)
            if (numel(handles.thresLevel)==4)
                separateThresholds                          = 1;
            else
                redGreenThresholds = input('Do you want to have different threshold levels for Green/Red Channels? Y/N [n]: ','s');
                if (strcmp(redGreenThresholds,'y'))||(strcmp(redGreenThresholds,'Y'))
                    separateThresholds                      = 1;
                else
                    separateThresholds                      = 0;
                end
            end
        else
            separateThresholds                              = 0;

        end
        if separateThresholds==1
            if numel(handles.minBlob)==1 
                handles.minBlob(2)=handles.minBlob(1); 
            end
            lowThresChangeG  = input (strcat('Low threshold  GREEN = ',num2str(handles.thresLevel(1)),'; '));
            highThresChangeG = input (strcat('High threshold GREEN = ',num2str(handles.thresLevel(2)),'; '));
            if numel(handles.thresLevel)==2
                handles.thresLevel(3) = handles.thresLevel(1);
                handles.thresLevel(4) = handles.thresLevel(2);
            end
            lowThresChangeR  = input (strcat('Low threshold  RED = ',num2str(handles.thresLevel(3)),'; '));
            highThresChangeR = input (strcat('High threshold RED = ',num2str(handles.thresLevel(4)),'; '));


            if ~isempty( lowThresChangeG);       handles.thresLevel(1) = lowThresChangeG; end
            if ~isempty( highThresChangeG);      handles.thresLevel(2) = highThresChangeG; end
            if ~isempty( lowThresChangeR);       handles.thresLevel(3) = lowThresChangeR; end
            if ~isempty( highThresChangeR);      handles.thresLevel(4) = highThresChangeR; end
            minBlob         = input (strcat('Minimum Size acceptable GREEN = [',num2str(handles.minBlob),']; '));
            if isempty( minBlob);      minBlob(1) = handles.minBlob(1); end
            
            minBlob2        = input (strcat('Minimum Size acceptable RED   = [',num2str(handles.minBlob),']; '));
            if isempty( minBlob2);      
                minBlob(2) = handles.minBlob(2);
            else
                minBlob(2) = minBlob2;
            end



        else
            lowThresChange  = input (strcat('Low threshold = ',num2str(handles.thresLevel(1)),'; '));
            highThresChange = input (strcat('High threshold = ',num2str(handles.thresLevel(2)),'; '));
            if ~isempty( lowThresChange);       handles.thresLevel(1) = lowThresChange; end
            if ~isempty( highThresChange);      handles.thresLevel(2) = highThresChange; end
            minBlob         = input (strcat('Minimum Size of Neutrophil acceptable = [',num2str(handles.minBlob),']; '));

        end
        rowsToPlot1     = input (strcat('Initial row to plot = [',num2str(indexR(1)),']; '));
        rowsToPlot2     = input (strcat('Final row to plot = [',num2str(indexR(end)),']; '));
        colsToPlot1     = input (strcat('Initial column to plot = [',num2str(indexC(1)),']; '));
        colsToPlot2     = input (strcat('Final column to plot = [',num2str(indexC(end)),']; '));

        if ~isempty( rowsToPlot1);
            initialRow = rowsToPlot1;
        else
            initialRow = indexR(1);
        end
        if ~isempty( rowsToPlot2);
            finalRow  = rowsToPlot2;
        else
            finalRow = indexR(end);
        end
        if ~isempty( colsToPlot1)
            initialColumn = colsToPlot1;
        else
            initialColumn = indexC(1);
        end
        if ~isempty( colsToPlot2);
            finalColumn = colsToPlot2;
        else
            finalColumn = indexC(end);
        end
        if ~isempty( minBlob);
            handles.minBlob = minBlob;
        end

        indexR = initialRow:finalRow;
        indexC = initialColumn:finalColumn;

        [handles] = verifyThresholdNeutrophils(handles,dataIn,indexR,indexC,indexL);
        close all
    else
        close all
    end
else
    close all;
end
