function firstNetwork=positionNeutrophil(dataIn,currFrameExternal,calculateSphericity)
%function firstNetwork=positionNeutrophil(dataIn,currFrameExternal,calculateSphericity)
%
%--------------------------------------------------------------------------
% positionNeutrophil  calculates the first metrics that will form handles.nodeNetwork:
%     first elements of the tracking process, get X, Y, Z for each neutrophil, then:
%       4      - distance to closest neighbour
%       5      - frame
%       10     - area of neutrophil
%       11     - label at frame 
%       15:20  - bounding box of neutrophil
%       other fields will be completed later on (6-ID ,7-parent, 8-child)
%
%       INPUT
%         dataIn:               path to labelled data
%
%         currFrameExternal:    currentFrame
%
%         calculateSphericity:  1 to calculate sphericity
%
%       OUTPUT
%         firstNetwork:         first network matrix as described above
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

firstNetwork=[];


if isa(dataIn,'char')
    % --------------------------- dataInName is not mat file, should be a folder with a)matlab files
    dir1                                    = dir(strcat(dataIn,'/*.mat'));
    numFrames                               = size(dir1,1);
    if ~exist('calculateSphericity','var')
        calculateSphericity=1;
    end
    for counterDir=1:numFrames
        tempDir                             = dir1(counterDir).name;
        dataInName                          =  strcat(dataIn,'/',tempDir);
        
        dataFromFile                        =  load(dataInName);
        
        if isfield(dataFromFile,'dataIn')
            dataIn2                         = dataFromFile.dataIn;
        elseif isfield(dataFromFile,'dataL')
            dataIn2 = dataFromFile.dataL;
        else
            namesF                          = fieldnames(dataFromFile);
            dataIn2                         = dataFromFile.dataL;
        end

        currentCentroids3                   = positionNeutrophil(dataIn2,counterDir,calculateSphericity);
        firstNetwork                        = [firstNetwork;currentCentroids3];  %#ok<AGROW>
    end
    if (size(firstNetwork,1) == 0)
       firstNetwork = zeros(0,28); 
    end
else
    numFrames=size(dataIn,4);
    for countFrames=1:numFrames
        %---- get properties of the labelled region
        currentProps                        =  regionprops(dataIn(:,:,:,countFrames),'centroid','area','boundingbox');
        
        if ~isempty(currentProps)
            [rows,cols,levs]                = size(dataIn);
            %-----  arrange centroid and area as a vertical vector
            currentCentroids                = [currentProps.Centroid]';
            currentArea                 	= [currentProps.Area]';
            if levs>1
                currentBox                      = reshape([currentProps.BoundingBox],[6 size(currentProps,1)])';
                %-----  centroids in [Y X Z]
                currentCentroids3               = [currentCentroids(2:3:end) currentCentroids(1:3:end) currentCentroids(3:3:end)];
            else
                currentBox                      = reshape([currentProps.BoundingBox],[4 size(currentProps,1)])';
                %-----  centroids in [Y X Z]
                currentCentroids3               = [currentCentroids(2:2:end) currentCentroids(1:2:end) ones(size(currentProps,1),1)];
            end
            %-----  distance between every element and all the rest
            [miuDist,dBPoints,numNeigh,compDist]     = distanceElements(currentCentroids3); %#ok<NASGU>
            
            %-----  compare with the minimum circular radius
            currentCentroids3(:,4)          =  miuDist;
            
            if ~exist('calculateSphericity','var')
                calculateSphericity=1;
            end
            if (exist('currFrameExternal','var'))&&(~isempty(currFrameExternal))
               currentCentroids3(:,5)       =  currFrameExternal;           
            else
                currentCentroids3(:,5)      =  countFrames;
            end
            currentCentroids3(:,10)         =  currentArea;
            %----- keep unique label at the current frame
            currentCentroids3(:,11)         = (1:size(currentArea,1));
            
            if levs>1
                currentCentroids3(:,15:20)      = currentBox;
            else
                currentCentroids3(:,[15 16 18 19])      = currentBox;
                currentCentroids3(:,17)      = 1;
                currentCentroids3(:,20)      = 1;
            end

            
            
            %----- a number to register how many neighbours at every 10 pix distance brackets
            currentCentroids3(:,28)         = numNeigh;
            %     %Obtain Sphericity of each neutrophil
            numNeutrop                                      = max(dataIn(:));
            if calculateSphericity==1
                for counterN = 1: numNeutrop
                    tempdata1                               = (dataIn==counterN);
                    if sum(tempdata1(:))>0
                        [rows,cols,levs]                                            = size(tempdata1);
                        indTemp                                 = find(tempdata1);
                        [XX,YY,ZZ]                              = ind2sub([rows cols levs],indTemp); %#ok<NASGU>

                        minR                                    = min(XX(:));
                        maxR                                    = max(XX(:));
                        minC                                    = min(YY(:));
                        maxC                                    = max(YY(:));

                        rowsCoor1                               = max(1,minR-1):min(rows,maxR+1);
                        colsCoor1                               = max(1,minC-1):min(cols,maxC+1);
                        
                        %calculate the surface by the zero crossings per slice
                        %surfNeutrop                             = zerocross(double(tempdata1(rowsCoor1,colsCoor1,:))-0.5);
                        %totSurfaceNeutrop                       = sum(surfNeutrop(:));
                        
                        %calculate the surface by finding the isosurface of the current object
                        % place in a container (tempdata2)larger than the object to avoid errors in the borders 
                        tempdata2                           = zeros(numel(rowsCoor1)+2,numel(colsCoor1)+2,levs+2);
                        tempdata2(2:end-1,2:end-1,2:end-1)  = tempdata1(rowsCoor1,colsCoor1,:);
                        %surfNeutrop                             = bwdist(1-tempdata1(rowsCoor1,colsCoor1,:));
                        %totSurfaceNeutrop2                       = sum(surfNeutrop(:)==1);
                        fv                                  = isosurface(tempdata2,0.5);
                        
                        totSurfaceNeutrop                   = size(fv.vertices,1);
                        volumeNeutrop                       = size(find(tempdata1(rowsCoor1,colsCoor1,:)),1);

                        
                        volSurfRatio                            = volumeNeutrop/totSurfaceNeutrop;
                        Sphericity                              = ((36*pi*(volumeNeutrop^2))^(1/3))/totSurfaceNeutrop;
%
                        currentCentroids3(counterN,26)          = volSurfRatio;
                        currentCentroids3(counterN,27)          = Sphericity;
                    else
                        currentCentroids3(counterN,26)          = -1;
                        currentCentroids3(counterN,27)          = -1;

                    end

                end

            end

            %----- append if applicable
            firstNetwork                = [firstNetwork;currentCentroids3];  %#ok<AGROW>
        end
        clear currentProps currentCe* currentA*;
    end
end
end


