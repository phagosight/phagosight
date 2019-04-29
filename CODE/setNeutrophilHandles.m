function handles = setNeutrophilHandles(dataIn,handles,thresLevels,askUser)
%function handles = setNeutrophilHandles(dataIn,handles)
%function handles = setNeutrophilHandles(dataIn,handles,thresLevels)
%function handles = setNeutrophilHandles(dataIn,handles,[],askUser)
%
%
%--------------------------------------------------------------------------
% setNeutrophilHandles  determines some characteristics of the
% datasets and assigns new fields
%                       to handles: 
%                               handles.rows,
%                               handles.cols,
%                               handles.levs,
%                               handles.numFrames
%                               handles.minBlob
%                               handles.thresLevel
%
%       INPUT
%         dataIn:           path to folder containing data or 3-D data matrix
%         handles:          handles struct
%         thresLevels:      1x2 vector containing minimum and maximum
%                           values to threshold data.
%         askUser:          1 for ask user to accept calculated threshold
%                           values
%
%       OUTPUT
%         handles:          handles struct updated as described above
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


if isa(dataIn,'char')
    % a name where the data is stored
    dir1 = dir(strcat(dataIn,'/*.mat'));
    tempDir = dir1(1).name;
    dataInName = strcat(dataIn,'/',tempDir);
    dataFromFile = load(dataInName);
    if isfield(dataFromFile,'dataIn')
        dataIn2 = dataFromFile.dataIn;
    else
        namesF = fieldnames(dataFromFile);
        dataIn2 = getfield(dataFromFile,namesF{1});
    end
    handles.numFrames = size(dir1,1); 
    [handles.rows,handles.cols,handles.levs] = size(dataIn2);
    if (strfind(dataIn,'_La'))
        % When the data input is _mat_La, there is no need to
        % calculate thresholds or data distributions, assume
        % everything is green and pass to the tracking.
        handles.ChannelDistribution=[1 handles.levs 0 0 0 0]';
        return;
    end    
    
else
    %a matlab matrix
    dataIn2=dataIn;
    if ~isfield(handles,'rows')
       [handles.rows,handles.cols,handles.levs,handles.numFrames]=size(dataIn2);
    end
end

if ~isfield(handles,'minBlob')
    handles.minBlob = 60;
end
    minData = min(dataIn2(:));
    maxData = max(dataIn2(:));
if (exist('thresLevels','var'))&&(numel(thresLevels)==0)
    disp('Manual Levels is empty')
    clear thresLevels;
end
if (exist('thresLevels','var'))&&((min(thresLevels)<minData)||...
                                  (max(thresLevels)>maxData))
    disp('Manual Levels do not match the levels of the input data')
    clear thresLevels;
end

if (exist('thresLevels','var'))
    % Manual levels that override Otsu
    if numel(thresLevels) == 1
        handles.thresLevel = [thresLevels 1.05*(thresLevels)];
        
    elseif numel(thresLevels)==2
        handles.thresLevel = [min(thresLevels) max(thresLevels)];
    else
        handles.thresLevel = ...
            [min(thresLevels(1:2)) max(thresLevels(1:2)) ...
             min(thresLevels(3:4)) max(thresLevels(3:4)) ];
    end
else
    %% get a dataR to calculate the threshold levels
    if isa(dataIn,'char')
        if (strfind(dataIn,'_Re'))
            dataR = double(dataIn2);
        else
            %try to read the _mat_Re dir to get the reduced version
            %of the data, if not, then reduce
            try
                dataRe = strcat(dataIn(1:end-2),'Re');
                dir2 = dir(strcat(dataRe,'/*.mat'));
                tempDir2 = dir2(1).name;
                dataInName2 = strcat(dataRe,'/',tempDir2);
                dataFromFile2 = load(dataInName2);
                dataR = dataFromFile2.dataR;
            catch
                dataR = reduceu(double(dataIn2),2);
            end
        end
    else
        dataR = double(dataIn);
    end
    
    %% Restrict data to fluorescence
    % If the Channels have been analysed and
    % handles.ChannelDistribution exists, then only use this to calculate the
    % levels. In the future, refine to obtain separate thresholds for green / red
    if isfield(handles,'ChannelDistribution')&&...
            (~isempty(handles.ChannelDistribution))
        indexLevels = ...
            [handles.ChannelDistribution(1):handles.ChannelDistribution(2) ...
                     handles.ChannelDistribution(5):handles.ChannelDistribution(6)];
        indexLevels (indexLevels==0) = [];
        dataR = dataR(:,:, indexLevels);
    else
        return;
    end
    
    %% Otsu's levels for thresholding
    dataIn3 = (double(dataR(:)));
    minData = (min(dataIn3));
    maxData = max((dataIn3)-minData);
    modeData = mode(dataIn3(:));
    meanData = mean(dataIn3(:));
    stdData = std(dataIn3(:));
    %%
    tempThres3(numel(indexLevels)) = 0;
    tempThres4(numel(indexLevels)) = 0;
    
    for counterL = 1:numel(indexLevels)
        %process the threshold at each image
        % take one slice at a time and filter to remove noise and
        % define better the thresholds
        tempSlice = (dataR(:,:,(counterL)));
        tempSlice2 = imfilter(tempSlice,gaussF(7,7,1),'replicate');
        % remove minimum and divide by maximum so that the
        % distribution is between [0 1]
        tempSlice3 = (tempSlice2-min(tempSlice2(:)));
        tempSlice3 = tempSlice3/max(tempSlice3(:));
        % Use OTSU to obtain threshold
        tempThres2 = graythresh (tempSlice3);
        if tempThres2==0
            %the current slice had no data
            tempThres3(counterL) = 0;
            tempThres4(counterL) = 0;
        else
            %to recover the thresholds from the data, binarise the
            %[0 1] data and use that to recover original
            %intensities, then discard zeros
            tempSlice4 = ((tempSlice3>tempThres2).*(dataR(:,:,counterL)));
            tempSlice5 = tempSlice4(:);
            tempSlice5(tempSlice5==0) = [];
            % get mean and std of those values that were above OTSU's level
            tempThres3(counterL) = mean(tempSlice5);
            tempThres4(counterL) = std(tempSlice5);          
        end
    end
    %% section added for ISBI Challenge to determine the levels and minBlob 
    if numel(tempThres3) ==1
        %there was just one slice, the higher threshold will be mean+std
        tempThres (2)  = tempThres3+tempThres4;
        tempThres3 (2) = tempThres3+tempThres4;
    else
        tempThres3 = sort (tempThres3);
        tempThres4 = sort (tempThres4);    
        tempThres = [tempThres3(end) tempThres3(end)+tempThres4(end)];
    end

    A           = tempThres3(end-1:end);% [level2 level1]
    alpha_G     = 0.180608;
    beta_G      = 0.251385;

    lt = modeData + alpha_G*(A(1) - minData);
    ht = modeData + beta_G*(A(2) - minData);
    handles.thresLevel = [lt ht];
    
    % Up to here is a different way of calculating levels, now calculate minimum
    % acceptable elements
    
    % Default values
    if handles.levs > 1
        % We have a 3D dataset.
        minimB = 110;
    else
        % for a 2D Data set
        minimB = 60;
    end
    % do a pre-segmentation if only one channel
    try
        if handles.ChannelDistribution(6)==0
            GreenLowT                       = bwlabeln(dataR(:,:,indexLevels)>lt);
            neutropToKeep                   = unique(GreenLowT.*(dataR(:,:,indexLevels)>ht));
            [dataL_green,numNeutropGreen]   = bwlabeln(ismember(GreenLowT,neutropToKeep(2:end)));
            
            regS = regionprops(dataL_green);
            
            M = mean([regS.Area]);
            S = std([regS.Area]);
            
            mB = max(minimB,(M-3*S)/4);
            
            handles.minBlob = mB;
        else
            handles.minBlob = minimB;
        end
    catch
        handles.minBlob = minimB;
    end
    %%
    %for very noisy data, the thresholds are too low, this will be indicated by
    %more than 30% of voxels above the threshold, then raise thresholds to minimum of
    %the mean + std,
    if (mean(dataIn3>handles.thresLevel(2))>0.3)
        handles.thresLevel = tempThres(end-1:end);% [level2 level1]
    end
    %it is necessary to pass dataIn in case there is red/green so
    %that they are properly identified by the indices
    
    %Modified to manage the variable askUser in verifyThresholdNeutrophils
    if(~exist('askUser','var'))
        askUser = true; 
    end
    if (askUser)
        [handles] = verifyThresholdNeutrophils(handles,dataIn); 
    end
    
end
