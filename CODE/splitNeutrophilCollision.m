function [numObjectsToSplit,dataL,handles] = splitNeutrophilCollision(handles,dataRe)
%function [numObjectsToSplit,dataL,handles] = splitNeutrophilCollision(handles,dataRe)
%
%--------------------------------------------------------------------------
% splitNeutrophilCollision  a routine to split all those cells that have been detected as a 
%     collision, will have to loop all the dataRs the results are SAVED to the 
%     files in the pre-existing folders 
%
%       INPUT
%         handles:              a struct with all parameters 
%         dataRe:               a string with the path to folder with dataR
%
%       OUTPUT
%         numObjectsToSplit:    number of split nodes
%         dataL:                path to labelled data folder
%         handles:              updated handles struct
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

% 
%%
% Generate the FILE NAME for the labelled folders
dataLa                                           = strcat(dataRe(1:end-2),'La');
%% Find all the objects that require splitting, keep them arranged according to the number of merged objects

% Generate the collision network based on the handles of the tracked data
collisionNetwork                                = detectNeutrophilCollision(handles);
% from the collision network determine which objects require splitting
collisionInTracks                               = collisionNetwork>1;
numObjectsToSplit                               = sum(collisionInTracks(:));
if numObjectsToSplit>0
    objectsToSplit                                  = zeros(numObjectsToSplit,2);
    indexObjects                                    = (find(collisionNetwork>1));
    % and the tracks to which the objects belong 
    
    %Objects to split (1-num object from nodeNetwork ,2-how many merged objects, 3-in which frame it exists)
    objectsToSplit(:,1)                             = (handles.finalNetwork(indexObjects));
    objectsToSplit(:,2)                             = (collisionNetwork(indexObjects));
    [objectsToSplit(:,3),indexObjsSorted]           = sort(objectsToSplit(:,1));
    objectsToSplit(:,4)                             = objectsToSplit(indexObjsSorted,2);
    objectsToSplit(:,1:2)                           = [];
    objectsToSplit(:,3)                             = handles.nodeNetwork(objectsToSplit(:,1),5);
    objectsToSplit(:,4)                             = handles.nodeNetwork(objectsToSplit(:,1),13);
    
    %% Find the frames where the objects exist
    FramesWithObjects                               = unique(objectsToSplit(:,3));
    FramesWithObjects(:,2)                          = histc(objectsToSplit(:,3),FramesWithObjects);
    
    numFramesToProcess                              = size(FramesWithObjects,1);
    
    %% Loops:  over the **FRAMES**  or loop over the **TRACKS**
    % Looping over the frames has the advantage that the dataR and dataLa are loaded once and it is not
    % necessary to re load them in other cases. BUT that does not make use of the knowledge of the
    % previous/posterior position of the cells when they are split. Looping over the tracks will require to
    % load dataR and dataLa several times but it will go with the cells and split them with knowledge of the
    % history of the cells. It could be that ALL frames that will be processed are read into an array, but
    % that most surely will have memory problems with larger data sets.
    
    % --------------------------- dataReName is not mat file, should be a folder with a)matlab files
    dir1                                            = dir(strcat(dataRe,'/*.mat'));
    dir2                                            = dir(strcat(dataLa,'/*.mat'));
    
    filtG                                           = gaussF(3,3,1);
    %% Loop over the FRAMES
    for counterDir=1:numFramesToProcess
        %%
        % This is the frame that requires processing
        tempDirRe                                   = dir1(FramesWithObjects(counterDir,1)).name;
        tempDirLa                                   = dir2(FramesWithObjects(counterDir,1)).name;
        dataReName                                  =  strcat(dataRe,'/',tempDirRe);
        dataLaName                                  =  strcat(dataLa,'/',tempDirLa);
        
        % Read first the intensity data "dataR"
        dataFromFile                                =  load(dataReName);
        if isfield(dataFromFile,'dataR')
            dataRe2                                 = dataFromFile.dataR;
        else
            namesF                                  = fieldnames(dataFromFile);
            dataRe2                                 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
        end
        % now read the labelled data dataL
        dataFromFile                                =  load(dataLaName);
        if isfield(dataFromFile,'dataL')
            dataLa2                                 = dataFromFile.dataL;
        else
            namesF                                  = fieldnames(dataFromFile);
            dataLa2                                 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
        end
        
        % dataL contains the previous segmentation, some of these objects need to be re-split according to
        % the intensity contained in dataR. Handles has the bounding box to process only a small area
        [rows,cols,levs]                            = size(dataRe2);
        [rowsL,colsL,levsL]                         = size(dataLa2);
        currObjects                                 = objectsToSplit(objectsToSplit(:,3)==FramesWithObjects(counterDir,1),1);
        
        %% smooth the current intensity
        dataRe3                                     = zeros(size(dataRe2));
        for counterLevs= 1:(levs/2)
            dataRe3(:,:,counterLevs)=imfilter(dataRe2(:,:,counterLevs),filtG);
        end
        %%
        % loop over the objects to split in the present Frame
        for k=1:FramesWithObjects(counterDir,2)
            %%
            labelObjectToSplit                      = handles.nodeNetwork(currObjects(k),11);
            %Proceed to split the labelled data
            initC                                   = max(1,-1+floor(handles.nodeNetwork(currObjects(k),15)));
            finC                                    = min(cols,1+ceil(handles.nodeNetwork(currObjects(k),15)+handles.nodeNetwork(currObjects(k),18)));
            initR                                   = max(1,-1+floor(handles.nodeNetwork(currObjects(k),16)));
            finR                                    = min(rows,1+ceil(handles.nodeNetwork(currObjects(k),16)+handles.nodeNetwork(currObjects(k),19)));
            currIntensity                           = ((dataRe3(initR:finR,initC:finC,1:levsL)));
            
            currLabelledFrame                       = ((dataLa2(initR:finR,initC:finC,1:levsL)));
           
            currLabelled                            = (currLabelledFrame==labelObjectToSplit);
                        
            currIntensLab                           = currIntensity.*(currLabelled);
            %%
            newSplit                                = watershed(-(currIntensLab));
            numDesiredObjects                       = objectsToSplit(objectsToSplit(:,1)==currObjects(k),2);
            
            if max(newSplit(:))==1
                % In the case that the partition was not achieved (only one region detected) try using the data without smoothing,
                % and then re-applying the watershed, if that fails, split in two randomly    
                currIntensity                           = ((dataRe2(initR:finR,initC:finC,1:levsL)));
                
                currIntensLab                           = currIntensity.*(currLabelled);
                
                %%
                newSplit                                = watershed(-(currIntensLab));
            end
            highestCurrentLabel                     = max(dataLa2(:));
            if max(newSplit(:))==1
                cumulObjects1=(cumsum(sum(sum(currLabelled,3),1))/sum(currLabelled(:)));
                cumulObjects2=(cumsum(sum(sum(currLabelled,3),2))/sum(currLabelled(:)));
                switch numDesiredObjects
                    case 2
                        %passing a 2 to all those on the right hand side
                        indexSide = find(cumulObjects1>0.5);
                        currLabelledFrame(:, indexSide,:)    = (currLabelledFrame(:,indexSide ,:)>0)*(highestCurrentLabel +1);
                        %
                    case 3
                        %partition into 3 horizontal sections
                        indexSide = find(cumulObjects1>0.33);
                        currLabelledFrame(:,indexSide ,:)    = (currLabelledFrame(:,indexSide ,:)>0)*(highestCurrentLabel +1);
                        indexSide = find(cumulObjects1>0.66);
                        currLabelledFrame(:,indexSide ,:)    = (currLabelledFrame(:,indexSide ,:)>0)*(highestCurrentLabel +2);
                        %
                    otherwise
                        %partition in 4 cuadrants
                        indexSideR = find(cumulObjects2<=0.5);
                        indexSideC = find(cumulObjects1>0.5);
                        currLabelledFrame(indexSideR,indexSideC ,:)      = (currLabelledFrame(indexSideR,indexSideC ,:)>0)*(highestCurrentLabel +1);
                        
                        indexSideR = find(cumulObjects2<=0.5);
                        indexSideC = find(cumulObjects1<=0.5);
                        currLabelledFrame(indexSideR,indexSideC  ,:)     = (currLabelledFrame(indexSideR,indexSideC ,:)>0)*(highestCurrentLabel +2);
                        
                        indexSideR = find(cumulObjects2>0.5);
                        indexSideC = find(cumulObjects1>0.5);
                        currLabelledFrame(indexSideR,indexSideC  ,:)      = (currLabelledFrame(indexSideR,indexSideC ,:)>0)*(highestCurrentLabel +3);
                        %
                end
                                
            else
                try
                    newLabelled                             = (double(newSplit).*currLabelled)>0;
                catch
                    q=1;
                end
                %% hierarchical Rejoining of objects
                % Once the current object has been split into MANY disjoint objects with a watershed transform, join
                % iteratively until just the desired number of objects is reached

                [disjointObjs,numDisjointObjects]       = bwlabeln(newLabelled);
                [disjointObjsProp]                      = regionprops(disjointObjs,'Area','Centroid');
                boundariesBetDisjoints                  = (disjointObjs==0).*currLabelled;


                % Find the distances between the centroids of the disjoint objects and use to link with COMPLETE linking
                clear X Y Z T
                X(numDesiredObjects,3)                  = 0; %#ok<AGROW>
                for k2 = 1:numDisjointObjects
                    if levsL>1
                        X(k2,:)                         = disjointObjsProp(k2).Centroid; %#ok<AGROW>
                    else
                        X(k2,1:2)                       = disjointObjsProp(k2).Centroid; %#ok<AGROW>
                        X(k2,3)                         = 1; %#ok<AGROW>

                    end
                end
                Y                                       = pdist(X);
                Z                                       = linkage(Y,'complete');
                % collapse to the number of objects that are required from this segmentation
                T                                       = cluster(Z,'maxclust',numDesiredObjects);

                % now assign the boundaries to the objects previously defined, sequentially first gets lucky
                finalLabel                              = zeros(size(disjointObjs));
                for k3=1:numDesiredObjects
                    currCluster                         = ismember(disjointObjs,find(T==k3));
                    
                    currBoundary                        = boundariesBetDisjoints.*imdilate(currCluster,ones(3,3,3));
                    
                    finalLabel                          = finalLabel+ k3*(currCluster + currBoundary);
                    boundariesBetDisjoints              = boundariesBetDisjoints&(1-currBoundary);
                end

                %finally, change the original label to a new label according to the new segmentation
                % what is 1 in finalLabel will remain with its current label, modify 2,3,...

                for k4=2:numDesiredObjects
                    currLabelledFrame(finalLabel==k4)    = highestCurrentLabel +k4 -1;
                end
            end
            dataLa2(initR:finR,initC:finC,1:levsL)      = currLabelledFrame ;
        end
        %%
        dataL = dataLa2; 
        save(dataLaName,'dataL');
        
    end
end

