function handles = trackNeutrophil(dataIn,handles,calculateSphericity)
%function handles = trackNeutrophil(dataIn,handles,calculateSphericity)
%
%--------------------------------------------------------------------------
% trackNeutrophil   receives the data as segmented neutrophils and then creates 
%     the concatenation of cell to cell between time frames, in the process, it 
%     calculates many other things like the velocity from a 3D sequence of frames 
%     for Neutrophils
%
%       INPUT 
%         dataIn:               path to folder with labelled data (3D frames
%                               stored in files)
%
%         handles:              handles struct
%
%         calculateSphericity:  1 to calculate the sphericity
%
%       OUTPUT
%         handles:              updated handles struct that include the field finalNetwork 
%                               [ X  Y  Z  distanceToNext timeSlice ID Parent Child ]
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



    %------ no input data is received, error -------------------------
    if nargin < 1
        help trackNeutrophil;
        handles=[]; return;
    end;
    %------ extra variable, to be used only when the tracking is done a second time, this allows
    %------ the extra linking of the tracks and deletion of smaller tracks -------------------
     if ~exist('Relink','var')
         ReLink = 0;
     end

    if ~exist('calculateSphericity','var')
        calculateSphericity=1;
    end
    %%

    %----- first, for every LABELLED  slice in time, Obtain centroids, area and eccentricity and return in a proper format
    %----- [1-xpos 2-ypos 3-zpos 4-miuDist 5-Frame 6-sequential 7-parent 8-child 9- 10-area 11-label at frame 12-whichRegion 13- track 14-Final Label 15-20 Bounding box ]
    firstNetwork                            = positionNeutrophil(dataIn,[],calculateSphericity);

    %----- define a circular radius for the cases when there is no parent and there is only one cell
    circRadius                              = ceil( (min(handles.rows,handles.cols))/4);
    % compare with the minimum circular radiuscurrentCentroids3(:,4)=
    % Unsure why this minimum comparison is performed keep the distances in a new Variable and replace at the
    % end:
    distanceToNeighbours                    = firstNetwork(:,4);
    firstNetwork(:,4)                       = min(firstNetwork(:,4)/2,circRadius);

    %----- Next assign parent / child relationships between the nodes
    if size(firstNetwork,1)>1
        %this will become handles.nodeNetwork
        parentNetwork                       = assignParentLoop(firstNetwork,handles);
        parentNetwork                       = assignParentOverlap(parentNetwork);
        %----- Reassign the concatenated nodes as branches of a bigger network. 
        % this will become handles.finalNetwork   
        neighbourhoodNetwork                = concatenateNodes(parentNetwork);

        %----- it may be the case that two nodes have been assigned the same child (most likely a merger will
        %----- occur, only one will keep it as tracks are assigned, so erase child from node that does not
        %----- get it and indicate that the nodes are possible mergers in column 32
        %%
        numObjects                              =size(parentNetwork,1);

        for k=1:numObjects
            indexRep                            = find(parentNetwork(:,8)==k);
            if numel(indexRep)>1
                parentNetwork(indexRep,32)      = 1;
            end
        end
        parentNetwork(:,8)                      = 0;
        for k=numObjects:-1:1
            indexPar = parentNetwork(k,7);
            if indexPar~=0
                parentNetwork(indexPar,8)=k;
            end
        end

        %%
        %----- clean the network from first nodes that could be incorrect
        reducedNetwork                      = cleanNetwork(neighbourhoodNetwork,firstNetwork);

        %----- link tracks that may have been oversplit
        linkedNetwork = reducedNetwork;
        if ~(ReLink == 1)
            finalNetwork                       = linkedNetwork;
        end
        %----- obtain distances of branches, average, total and per hop
        if ~isempty(finalNetwork)
            distanceNetwork                 = getDistanceNet(finalNetwork,parentNetwork);

            %----- remove outliers i.e. those branches whose avDistance > mean(avDistance)+3*std(avDistance)
            mm                              = mean(distanceNetwork.avPerTrack);
            ss                              = std(distanceNetwork.avPerTrack);

            [q1,q2]                         = find(distanceNetwork.avPerTrack>(mm+3*ss));

            distanceNetwork.perHop(:,q2)    = [];
            distanceNetwork.numHops(q2)     = [];
            distanceNetwork.totPerTrack(q2) = [];
            distanceNetwork.maxPerTrack(q2) = [];
            distanceNetwork.avPerTrack(q2)  = [];
            finalNetwork(:,q2)              = [];

            distanceNetwork.biasedMean      = mean(distanceNetwork.avPerTrack);
            distanceNetwork.biasedStd       = std(distanceNetwork.avPerTrack);

            %----- check that neutrophils have common parents to determine overlap
            finalLabel                          = neutrophilOverlap(firstNetwork,finalNetwork,handles);

        else
            finalNetwork                    = [];
            distanceNetwork.perHop          = [];
            distanceNetwork.numHops         = [];
            distanceNetwork.totPerTrack     = [];
            distanceNetwork.maxPerTrack     = [];
            distanceNetwork.avPerTrack      = [];
            distanceNetwork.biasedMean      = 0;
            distanceNetwork.biasedStd       = 0;
            finalLabel                          =[];

        end



    else
        parentNetwork                       = [];
        neighbourhoodNetwork                = [];
        reducedNetwork                      = [];
        linkedNetwork                       = [];
        finalNetwork                        = [];
        distanceNetwork.perHop              = [];
        distanceNetwork.numHops             = [];
        distanceNetwork.totPerTrack         = [];
        distanceNetwork.maxPerTrack         = [];
        distanceNetwork.avPerTrack          = [];
        distanceNetwork.biasedMean          = 0;
        distanceNetwork.biasedStd           = 0;
    end
    % Assign all the network parameters to the handles

    parentNetwork(:,4)                      = distanceToNeighbours;
    handles.nodeNetwork                     = parentNetwork;
    handles.neighNetwork                    = neighbourhoodNetwork;
    handles.reducedNetwork                  = reducedNetwork;
    handles.linkedNetwork                   = linkedNetwork;
    handles.finalNetwork                    = finalNetwork;
    handles.distanceNetwork                 = distanceNetwork;
    %handles.finalNetwork(:,q2)=[];
    handles.finalLabel                      = finalLabel;

    %% assign to which track they belong (13), a final label (14), and velocity (9)
    for counterTracks=1:size(handles.finalNetwork,2)
        trackIndex                                  = find(handles.finalNetwork (:,counterTracks));
        nodeIndex                                   = handles.finalNetwork(trackIndex,counterTracks); %#ok<FNDSB>
        
        handles.nodeNetwork(nodeIndex,13)           = counterTracks; 
        handles.nodeNetwork(nodeIndex,14)           = handles.finalLabel(counterTracks); 
        handles.nodeNetwork(nodeIndex(2:end),9)     = sqrt(sum(diff(handles.nodeNetwork(nodeIndex,1:3)).^2,2));
    end

    %% Last part, assign distance to appearing/disappearing neighbours

    for counterF =1:handles.numFrames

        [miuD,distanceBetPoints,numNeighboursAtDist,compDistance,distToDisappearing,distToAppearing] =distanceElements (handles.nodeNetwork(handles.nodeNetwork(:,5)==counterF,[1 2 3 7 8]));

        %----- distance to those objects that Appear or Disappear at a certain frame if none, then
        %      distance is INF
        if isempty(distToDisappearing)
            distToDisappearing                      = inf;
        end
        if isempty(distToAppearing)
            distToAppearing                         = inf;
        end
        handles.nodeNetwork(handles.nodeNetwork(:,5)==counterF,30) = distToDisappearing;
        handles.nodeNetwork(handles.nodeNetwork(:,5)==counterF,31) = distToAppearing;

    end

