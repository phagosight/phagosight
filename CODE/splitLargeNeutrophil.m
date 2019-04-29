function [handles,numObjectsToSplit] = splitLargeNeutrophil(handles,dataRe,...
                                                  outlierThreshold)
%function [handles,numObjectsToSplit] = splitLargeNeutrophil(handles,dataRe,...
% outlierThreshold)
%
%--------------------------------------------------------------------------
% splitLargeNeutrophil  a routine to split all those cells that
% have been detected as a very large outlier, will have to loop all
% the dataRs the results are SAVED to the files in the pre-existing folders.
%
%       INPUT
%         handles:              a struct with all parameters 
%         dataRe:               a string with the path to folder with dataR
%         outlierThreshold:     the outlier threshold
%
%       OUTPUT
%   	  handles:              updated handles
%         numObjectsToSplit:    number of split object 
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
% method in the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

% a routine to split all those cells that have been detected as a
% collision, will have to loop all the dataRs
%%
%% Find all the objects that require splitting, keep them arranged
%% according to the number of merged objects

if ~isfield(handles,'nodeNetwork')
    % positionNeutrophil calculates sphericity if third argument is
    % set to 1
    firstNetwork = positionNeutrophil(dataRe,[],0);
    numNodesInNet = size(firstNetwork,1);
    if (numNodesInNet > 0)
        firstNetwork(:,6) = (1:numNodesInNet)';
        handles.nodeNetwork = firstNetwork;
    else
        handles.nodeNetwork = [];
        numObjectsToSplit = 0;
        return;
    end
end
% Generate the FILE NAME for the labelled folders
%%
if strcmp(dataRe(end-1:end),'Re')
    dataLa = strcat(dataRe(1:end-2),'La');
elseif strcmp(dataRe(end-1:end),'La')
    dataLa = dataRe;
    dataRe = strcat(dataRe(1:end-2),'Re');
end

if (~exist('outlierThreshold','var'))||(isempty(outlierThreshold))
%%

    indexGreen = (handles.ChannelDistribution(1):handles.ChannelDistribution(2));
    indexRed = (handles.ChannelDistribution(5):handles.ChannelDistribution(6));
    
    
    %%%%% Detection based on variable thresholds along time for the
    %%%%% time being same for red and green
    %calculate the volume at every time point, the STD and the     
    %mean + 4*STD but filter them 
    
    volAtTime = zeros(handles.numFrames,3);
    for counterTime = 1:handles.numFrames
        % Calculate the mean
        volAtTime(counterTime,1) = ...
            mean(firstNetwork(firstNetwork(:,5)==counterTime,10));        
        % Calculate the STD
        volAtTime(counterTime,2) = ...
            std(firstNetwork(firstNetwork(:,5)==counterTime,10));
    end
    %complete for boundaries
    volAtTime(volAtTime(:,2)==0,2) = 0.25*volAtTime(volAtTime(:,2)==0,1);
    
    %% IT IS ASSUMED THAT THERE IS A MINIMUM NUMBER OF FRAMES to
    %% avoid bias, but for small sets it will simplify the calculations
    if handles.numFrames>8
        %calculation of the outliers, use mean value for centre but
        %median for STD to avoid biasing
        for counterTime = 4:handles.numFrames-3
            volAtTime(counterTime,3) = ...
                mean(volAtTime(counterTime-3:counterTime+3,1))+...
                4*(median(volAtTime(counterTime-3:counterTime+3,2)))+1;
        end
        %complete the boundaries
        volAtTime(1:3,3) = max(volAtTime(4:5,3));
        volAtTime(end-2:end,3) = max(volAtTime(end-5:end-3,3));
        %extend the threshold per time to all points of network
        
        outlierThreshold = zeros(firstNetwork(end,6),1);
        for counterTime = 1:handles.numFrames
            outlierThreshold = outlierThreshold+...
                volAtTime(counterTime,3)*(firstNetwork(:,5)==counterTime);
        end
    else        
        outlierThreshold = mean(volAtTime(:,1))+4*(median(volAtTime(:,2)))+1;
    end
    objectsToSplit(:,1:3) = ...
        handles.nodeNetwork(handles.nodeNetwork(:,10)>outlierThreshold,[6 10 5]);
end

%Objects to split (1-num object from nodeNetwork ,2-how many merged
%objects, 3-in which frame it exists)
% Find the frames where the objects exist
FramesWithObjects = unique(objectsToSplit(:,3));
FramesWithObjects(:,2) = histc(objectsToSplit(:,3),FramesWithObjects);

numFramesToProcess = size(FramesWithObjects,1);
numObjectsToSplit = size(objectsToSplit,1);

% Loops:  over the **FRAMES**
% Looping over the frames has the advantage that the dataR and
% dataLa are loaded once and it is not necessary to re load them in
% other cases. BUT that does not make use of the knowledge of the
% previous/posterior position of the cells when they are
% split. Looping over the tracks will require to load dataR and
% dataLa several times but it will go with the cells and split them
% with knowledge of the history of the cells. It could be that ALL
% frames that will be processed are read into an array, but
% that most surely will have memory problems with larger data sets.

% ----------- dataReName is not mat file, should be a folder with
% a)matlab files

%disp('Split Very large cells and re-save in a the same folder _mat_La')
if (numFramesToProcess>1)
    disp('Split Large Neutrophils and re-save in same folder  *_mat_La')
end
dir1 = dir(strcat(dataRe,'/*.mat'));
dir2 = dir(strcat(dataLa,'/*.mat'));
filtG = gaussF(5,5,1);
%% Loop over the FRAMES
for counterDir=1:numFramesToProcess
    %%
    tempDirRe = dir1(FramesWithObjects(counterDir,1)).name;
    tempDirLa = dir2(FramesWithObjects(counterDir,1)).name;
    dataReName =  strcat(dataRe,'/',tempDirRe);
    dataLaName =  strcat(dataLa,'/',tempDirLa);
    
    % Read first the intensity data "dataR"
    dataFromFile =  load(dataReName);
    if isfield(dataFromFile,'dataR')
        dataRe2 = dataFromFile.dataR;
    else
        namesF = fieldnames(dataFromFile);
        dataRe2 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
    end
    % now read the labelled data dataL
    dataFromFile =  load(dataLaName);
    if isfield(dataFromFile,'dataL')
        dataLa2 = dataFromFile.dataL;
    else
        namesF = fieldnames(dataFromFile);
        dataLa2 = getfield(dataFromFile,namesF{1});
    end

    %process all those objects requiring to be split in a loop
    %
    [rows,cols,levs] = size(dataRe2);
    currObjects = ...
        objectsToSplit(objectsToSplit(:,3)==FramesWithObjects(counterDir,1),1);    
    dataRe3 = dataRe2;
    
    %%
    % loop over the objects to split in the present Frame
    for k=1:FramesWithObjects(counterDir,2)
        
        %Proceed to split the labelled data
        labelObjectToSplit = handles.nodeNetwork(currObjects(k),11);
        %detect the region where the merged cells occur, this
        %reduces time considerably
        %   Columns
        initC = max(1,-1+ floor(handles.nodeNetwork(currObjects(k),15)));
        finC = ...
            min(cols,1+ceil(handles.nodeNetwork(currObjects(k),15)+...
                            handles.nodeNetwork(currObjects(k),18)));
        %   Rows
        initR = max(1,-1+ floor(handles.nodeNetwork(currObjects(k),16)));
        finR = ...
            min(rows,1+ceil(handles.nodeNetwork(currObjects(k),16)+...
                            handles.nodeNetwork(currObjects(k),19)));
        % partition over the maximum intensity projection of EACH Current CHANNEL
        if max(handles.ChannelDistribution)==1
            %this is the case that there is only one channel
            initL=1;
            finL=1;            
        elseif ((handles.nodeNetwork(currObjects(k),3))>=...
                handles.ChannelDistribution(1))&&...
                ((handles.nodeNetwork(currObjects(k),3))<=...
                 handles.ChannelDistribution(2))
            %this belongs to the green channel
            initL = handles.ChannelDistribution(1);
            finL = handles.ChannelDistribution(2);
        else
            %this belongs to the red channel
            initL = handles.ChannelDistribution(5);
            finL = handles.ChannelDistribution(6);
        end
        
        currIntensity = dataRe3(initR:finR,initC:finC,initL:finL);
        currLabelledFrame = dataLa2(initR:finR,initC:finC,:);
        currLabelledFrame(currLabelledFrame~=labelObjectToSplit)=0;
        
        currIntensityT = (sum(currIntensity,3));
        
        currLabelled = (currLabelledFrame==labelObjectToSplit);
        
        currLabelledT = max(currLabelled,[],3);
        currIntensLabT = currIntensityT.*(currLabelledT);
        % find the watershed of the filtered data, BUT only use to
        % split if there is just one watershed and two clear regions
        newSplitT = watershed(-(imfilter(currIntensLabT,filtG)));
        
        numRegionsInSplit = max(newSplitT(:));
        highestCurrentLabel = max(dataLa2(:));
        %%
        if (numRegionsInSplit==2)
            % if there are still no boundaries then try by eroding the classes
            % numRegionsInSplit==1 OR there are too many, numRegionsInSplit>8
            
            clear ratioVol
            for counterRegion = 1:numRegionsInSplit
                % get areas of sides, then select the most balanced
                % split in terms of volume
                tempArea = (newSplitT==counterRegion);
                tempAreaC = bwlabel(newSplitT~=counterRegion,4);
                if max(tempAreaC(:))>1
                    %the area is central and the others appear at
                    %different sides, do not consider
                    ratioVol(counterRegion) = 0;
                else
                    %calculate the areas of the region and its complement
                    volArea = sum(sum(sum(...
                      currLabelled.*repmat(tempArea, [1 1 size(currLabelled,3)]))));
                    volAreaC = sum(sum(sum(...
                      currLabelled.*repmat(tempAreaC,[1 1 size(currLabelled,3)]))));
                    ratioVol(counterRegion) = ...
                        min(volArea,volAreaC)/max(volArea,volAreaC);
                end
            end
            [maxRatioVol,indMaxRatioVol] = max(ratioVol);
            regionToBoost = ...
                currLabelled.* repmat((newSplitT==indMaxRatioVol),...
                                      [1 1 size(currLabelled,3)]);
            regionToDelete = ...
                imdilate(regionToBoost,ones(3)).*repmat((newSplitT==0),...
                                                        [1 1 size(currLabelled,3)]);
            
            currLabelledFrame(regionToDelete>0) = 0;
            currLabelledFrame(regionToBoost>0) = highestCurrentLabel +1;
            %test if the remaining region has been split in between
            remainingRegions = bwlabeln(currLabelledFrame==labelObjectToSplit);
            numRemainingRegions = max(remainingRegions(:));
            if numRemainingRegions>1
                for k3=2:numRemainingRegions
                    currLabelledFrame(remainingRegions==k3) = ...
                        highestCurrentLabel +k3;                 
                end
            end
                    dataLa2(initR:finR,initC:finC,:) = currLabelledFrame ;
        else
            %eroding the classes and split into two classes
            distFromBackground = zeros(size(currLabelledFrame));
            VoxToDelete = zeros(size(currLabelledFrame));
            for counterL = 1:size(currLabelledFrame,3)
 
                tempLabelledSlice = (currLabelledFrame(:,:,counterL)>0);
                if max(tempLabelledSlice(:))>0
                    distFromBackground(:,:,counterL) = bwdist(1-tempLabelledSlice);
                end
            end
            %%            
            numClasses = 1;
            distToErode = 0;
            
            ratioAreaClasses = 0;
            try
                while (ratioAreaClasses<0.1)||(numClasses<2)
                    distToErode = distToErode+1;
                    [sepObjs, numClasses] = bwlabeln(distFromBackground>distToErode);
                    areaClasses = regionprops(sepObjs,'Area');
                    ratioAreaClasses = ...
                        min([areaClasses.Area])/max([areaClasses.Area]);
                end
                
                if distToErode<3
                    %%
                    distFromObject1 = bwdist(sepObjs==1);
                    distFromObject2 = bwdist(sepObjs==2);
                    VoxToDelete1 = ...
                        (distFromObject1<=(distToErode+2))+...
                        (distFromObject2<=(distToErode+2));
                    for counterL=initL:finL
                        VoxToDelete(:,:,counterL) = ...
                            imdilate(bwmorph(bwmorph(bwmorph(...
                                VoxToDelete1(:,:,counterL)==2,'skel',distToErode),...
                                                     'diag'),'bridge'),[1 1;1 1]);
                    end
                    
                    currLabelledFrameVal = unique(currLabelledFrame);
                    
                    finalObject1 = ...
                        max(currLabelledFrameVal)*...
                        ((distFromObject1<(distToErode))&(1-VoxToDelete));
                    finalObject2 = ...
                        (highestCurrentLabel +1)*...
                        ((distFromObject2<(distToErode))&(1-VoxToDelete));
                    currLabelledFrame = ...
                        finalObject1+(finalObject2.*(finalObject1==0));
                    %%                                       
                    dataLa2(initR:finR,initC:finC,:) = currLabelledFrame ;
                end
            catch
                %It was not possible to split through erosion
            end
           
        end
        
    end

    dataL = dataLa2;
    save(dataLaName,'dataL');

end
