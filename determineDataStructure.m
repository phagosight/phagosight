function [handles,dataRe2] = determineDataStructure(dataRe,dataStructure)
%function [handles,dataRe2] = determineDataStructure(dataRe)
%function [handles,dataRe2] = determineDataStructure(dataRe,dataStructure)
%--------------------------------------------------------------------------
% DetermineDataStructure    determines the distribution of channels from the 
%     data sets, and assigns them to the handles, it also provides a GUI to 
%     allow the user to verify the distribution
%
%       INPUT  
%         dataRe:           path to the reduced data.
%         dataStructure:    6x1 vector containing the slice's first and
%                           last indexes of the Green, DIC and Red Channel
%                           Example: [6 10 1 5 11 15]' means slices 1 to 5
%                           correspond to the DIC channel, slices 6 to 10
%                           to the Green channel and, finally, slices from
%                           11 to 15 to the Red channel.
%       OUTPUT  
%         handles:          handles structure containing the
%                           ChannelDistribution field.
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
if isa(dataRe,'char')
    %dataRe is the path where the data is stored
    dir1 = dir(strcat(dataRe,'/*.mat'));
    tempDir = dir1(1).name;
    dataReName = strcat(dataRe,'/',tempDir);
    dataFromFile = load(dataReName);
    if isfield(dataFromFile,'dataRe')
        dataRe2 = dataFromFile.dataRe;
    else
        namesF = fieldnames(dataFromFile);
        if isempty(namesF(contains(namesF,'dataR')))
            dataRe2 = dataFromFile.(namesF{1});
        else
            dataRe2 = dataFromFile.(namesF{contains(namesF,'dataR')});
        end
    end
    handles.numFrames = size(dir1,1); 
    [handles.rows,handles.cols,handles.levs] = size(dataRe2);
else
    %dataRe is a matlab matrix
    dataRe2=dataRe;
    [handles.rows,handles.cols,handles.levs,handles.numFrames]=size(dataRe2);
end


%% Basic data has been acquired, rows, cols, levs, frames, now determine how it is distributed

if (handles.levs==1)
    handles.ChannelDistribution             =  [1 1 0 0 0 0 ]';
    
elseif exist('dataStructure','var')
    
    handles.ChannelDistribution             = dataStructure(:);
else
    
    currHist(20,handles.levs)                                   = 0;
    maxData                                                     = double(max(dataRe2(:)));
    minData                                                     = double(min(dataRe2(:)));
    rangeData                                                   = round(linspace(minData,maxData,20));
    
    for counterLevs = 1:handles.levs
        currLevel                                               = dataRe2(:,:,counterLevs);
        currHist(:,counterLevs)                                  = hist(currLevel(:),rangeData);
        
        
    end
    
    
    %%
    [rows,cols,levs]= size (dataRe2);
    topPart = dataRe2(:,:,1:floor(levs/2));
    botPart = dataRe2(:,:,1+floor(levs/2):end);
    
    minT    = double(min(topPart(:)));
    maxT    = double(max(topPart(:)));
    minB    = double(min(botPart(:)));
    maxB    = double(max(botPart(:)));
    
    xAxis = linspace(min(minB,minT),max(maxT,maxB),100);
    
    [yT,xT] = hist(topPart(:),xAxis); %#ok<NASGU>
    [yB,xB] = hist(botPart(:),xAxis); %#ok<NASGU>

    
    similarityM         = sum(abs(yT-yB))/rows/cols/levs;
    similarityM2        = (sum(cumsum(yT)-cumsum(yB)));
    
    % pre assign a distribution thinking the data is distributed 50-50. presume also that the Red Channel is
    % not present
    if similarityM>0.5
        % high similarity measure implies that the top half contains the fluorescence and the bottom DIC
        if similarityM2>0
            % positive similarityM2 implies that the fluorescence in on the top
            greenInit   = 1;
            greenFin    = round(levs/2);
            DICInit     = round(levs/2)+1;
            DICFin      =       levs;
        else
            % DIC on top
            DICInit     = 1;
            DICFin      = round(levs/2);
            greenInit   = round(levs/2)+1;
            greenFin    =       levs;
        end
    else
        %Low similarity implies both sections have fluorescence
        greenInit   = 1;
        greenFin    = levs;
        DICInit     = 0;
        DICFin      = 0;
    end
    RedInit = 0;
    RedFin = 0;
    
    
    
    %% Generate the figure with confirmation of dimensions
    handleFigConf = figure(11);
    clf
    set(handleFigConf,'position',[233   289   1370   550 ],'menuBar','none','toolBar','none','numbertitle','off')
    h1=subplot(131);
    ribbon(currHist)
    set(h1,'ytick',([5 10 15 20]),'yticklabel',rangeData([5 10 15 20]),'position',[0.03766 0.42 0.26988 0.5]);
    view(110,40)
    rotate3d on;
    title('Histograms of individual images of one time frame','fontsize',14)
    xlabel('Levels')
    ylabel('Intensity level')
    zlabel('number of elements')
    
    textExplanation         = ' Please, confirm the distribution of the channels on the data levels (slices). The histogram above can be a guide as DIC concentrates in the centre and Fluorescence on the lower intensities. Type the correct initial and final slices of each channel and when ready press "OK".';
    textExplanationH        = uicontrol(handleFigConf,'Style','text','string',textExplanation  ,    'position',[48   130    330   100],'foregroundcolor',0.0*[1 1 1],'fontsize',10, 'tag','textDICChannel','backgroundcolor',0.8*[1 1 1]);
    set(textExplanationH,'position',[30 30 420 140],'horizontalAlignment','left','fontsize',13.5)
    %%
    h211=subplot(132);
    if (greenInit~=0)
        imagesc(dataRe2(2:end-1,2:end-1,greenInit))
        title('First slice of Green Channel','fontsize',14)
    else
        imagesc(repmat(dataRe2(2:end-1,2:end-1,1)/max(max(dataRe2(:,:,1))),[1 1 3]))
        title('First slice of z-Stack','fontsize',14)        
    end
    h212=subplot(133);
    if (DICInit~=0)
        imagesc(repmat(dataRe2(2:end-1,2:end-1,DICInit)/max(max(dataRe2(:,:,DICInit))),[1 1 3]))
        title('First slice of DIC Channel','fontsize',14)
    else
        imagesc(repmat(dataRe2(2:end-1,2:end-1,end)/max(max(dataRe2(:,:,end))),[1 1 3]))
        title('Last slice of z-Stack','fontsize',14)
    end
    %%
    qq=hot;
    qq2(:,2)=qq(:,1);
    qq2(:,3)=qq(:,3);
    qq2(:,1)=qq(:,3);
    colormap(qq2)

    
    %% Pass the values through an input dialog
    prompt                                  ={'Green Initial Slice','Green Final Slice','DIC Initial Slice ','DIC  Final Slice','Red Initial Slice ','Red Final Slice '};
    name                                    ='Enter Channel Levels';
    numlines                                =2*[1 1 1 1 1 1]';
    defaultanswer                           ={num2str(greenInit),num2str(greenFin),num2str(DICInit),num2str(DICFin),num2str(RedInit),num2str(RedFin)};
    
    options.Resize                          ='on';
    options.WindowStyle                     ='normal';
    options.Interpreter                     ='tex';
    answer                                  =inputdlg(prompt,name,numlines,defaultanswer,options);
    %%    
    handles.ChannelDistribution             =  str2double([answer]);
    close (handleFigConf);
end
