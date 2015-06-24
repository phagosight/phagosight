function parentNetwork = assignParentLoop(firstNetwork,handles)
%function parentNetwork = assignParentLoop(firstNetwork,handles)
%
%--------------------------------------------------------------------------
% assignParentLoop  routine that assigns parent-child relationship for the 
%     tracking process.
%
%     INPUT 
%         firstNetwork: a matrix that contains one row **per object** 
%             segmented. Columns correspond to:
%               *1,2,3      - X, Y, Z for each neutrophil
%               *4          - distance to closest neighbour
%               *5          - frame
%               *10         - area of neutrophil
%     OUTPUT
%         parentNetwork: a copy of the firstNetwork matrix containing 
%             new additional information about the parent-child 
%             relationship of the nodes in firstNetwork. Columns 6-ID,
%             7-parent and 8-child are completed here.
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
levs                                    = handles.levs;

%----- Dimension check, number of nodes in the network
[numNodesInNet]                         = size(firstNetwork,1);
%----- assign an ID for the nodes of the Network (some may be removed later on, so keep this to identify parent/child
firstNetwork(:,6)                       = (1:numNodesInNet)';
%----- use circular radius as a minimum case of distance for the other cases, when there is only one node per frame
meanDistance                            = mean(firstNetwork(:,4));
stdDistance                             = std (firstNetwork(:,4));
%----- Avoid very large circular radius, keep maximum  to    mu+std
firstNetwork(:,4)                       = min(floor(meanDistance+stdDistance),firstNetwork(:,4));

%----- Channel Distribution
indexGreen                                          = (handles.ChannelDistribution(1):handles.ChannelDistribution(2));
indexRed                                            = (handles.ChannelDistribution(5):handles.ChannelDistribution(6));


%----- keyHoleDef  will contain the X Y Z position of parent and child to be tested to generate KeyHole
for counNodes=1:numNodesInNet-1
    %----- start the assignation process by checking the area of the currrent object, 
    %----- This line will discard all nodes whose area is less than   X pixels AND those with
    %----- child
    if (firstNetwork(counNodes,10)>0)&&(firstNetwork(counNodes,8)==0)
        %----- select a time frame and process objects of this frame
        currentFrame                    = firstNetwork(counNodes,5);
        %----- has the node to be analysed a parent? It is important for the sequences,
        %-----     if YES, assign them for the tailCone (predict landing),
        %-----     if NO, repeat its position and a circle will applied
        if firstNetwork(counNodes,7)==0
            %----  no parent
            keyHoleDef                  = repmat(firstNetwork(counNodes,1:3),2,1);
        else
            %------ place parent BELOW child
            keyHoleDef                  = (firstNetwork(counNodes,1:3));
            k=counNodes;
            while firstNetwork(k,7)~=0
                % repeat for all parents of the sequence
                keyHoleDef              = ([keyHoleDef;firstNetwork(firstNetwork(k,7),1:3) ]);
                k                       = firstNetwork(k,7);
            end
        end
        
        %Determining if there were 2 fluorescent channels using indices of Green and Red
        if (max(indexGreen)>0)&&(max(indexRed)>0)
            %Determining parenhood and finding all objects that belong to 
            %the following frame to test for parenthood
            if (ismember(round(keyHoleDef(1,3)),indexGreen))
                %ParentGreen
                childLevel                      = firstNetwork((firstNetwork(:,5)==(1+currentFrame))&(ismember(round(firstNetwork(:,3)),indexGreen)),:);
            else
                %ParentRed
                childLevel                      = firstNetwork((firstNetwork(:,5)==(1+currentFrame))&(ismember(round(firstNetwork(:,3)),indexRed)),:);
            end
        else
            childLevel                      = firstNetwork(firstNetwork(:,5)==(1+currentFrame),:);
        end

        %------ Check track for possible children
        [inTheTail,distanceToParent,whichRegion] = testKeyHole(keyHoleDef(1:2,:),childLevel(:,1:3),firstNetwork(counNodes,4));
        %------ if there is any node in the tail assign relationships
        if any(inTheTail>0)
            %------- if there is more than one child node that falls in the footprint, then assign
            %------- paternity to the closest one, since the network is sequential, it should arrange them properly
            distTemp                    = distanceToParent;
            distTemp(inTheTail==0)      = inf;
            %------- this finds the closest RBC to the parent, should Area be used???
            [minD,indexD]               = min(distTemp);
            %------- if the index has been erased, the relationships have been assigned, skip
            if childLevel(indexD,6)~=0
                if (max(indexGreen)>0)&&(max(indexRed)>0)
                    if ismember(round((firstNetwork(counNodes,3))),indexGreen)
                        ParentBelow = 1;
                    else
                        ParentBelow = 0;
                    end
                    if ismember(round((firstNetwork(childLevel(indexD,6),3))),indexGreen)
                        ChildBelow = 1;
                    else
                        ChildBelow = 0;
                    end
                    %Only assign if both nodes are in the same half of the data
                    if ChildBelow==ParentBelow
                        %---- parenthood is assigned to CHILD NODE
                        firstNetwork(childLevel(indexD,6),7)    = firstNetwork(counNodes,6);
                        %---- childhood is assigned to  PARENT NODE
                        firstNetwork(counNodes,8)               = childLevel(indexD,6);
                        %---- region is assigned to  Parent NODE
                        firstNetwork(counNodes,12)               = whichRegion(indexD);
                    end
                else
                    %---- parenthood is assigned to CHILD NODE
                    firstNetwork(childLevel(indexD,6),7)    = firstNetwork(counNodes,6);
                    %---- childhood is assigned to  PARENT NODE
                    firstNetwork(counNodes,8)               = childLevel(indexD,6);
                    %---- region is assigned to  Parent NODE
                    firstNetwork(counNodes,12)               = whichRegion(indexD);
                end
                %---- erase index so the node cannot be reassigned
                childLevel(indexD,6)    = 0;
            end
        end
    end
end

parentNetwork=firstNetwork;
