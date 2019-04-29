function firstNetwork = assignParentOverlap(firstNetwork)
%function firstNetwork = assignParentOverlap(firstNetwork)
%
%--------------------------------------------------------------------------
% assignParentLoop  routine that assigns parent-child relationship (based 
%     on overlapping between nodes) for the tracking process.
%
%     INPUT 
%         firstNetwork: a matrix that contains one row **per object** 
%             segmented. Columns correspond to:
%               *1,2,3      - X, Y, Z for each neutrophil
%               *4          - distance to closest neighbour
%               *5          - frame
%               *10         - area of neutrophil
%     OUTPUT
%         firstNetwork: a copy of the firstNetwork matrix containing 
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

%----- Dimension check, number of nodes in the network
[numNodesInNet]                         = size(firstNetwork,1);
firstNetwork(:,6)                       = (1:numNodesInNet)';
for counNodes=1:numNodesInNet-1
    %This procedure will only be applied to nodes without children
    if firstNetwork(counNodes,8)==0
        %----- start the assignation process with the simplest overlap rule:
        %-----     LARGE overlap between bounding boxes AND neighbours far away
        currentFrame                       = firstNetwork(counNodes,5);
        nextFrame                          = currentFrame+1;
        childLevel                         = firstNetwork((firstNetwork(:,5)==(nextFrame)),:);
        numChild                           = size(childLevel,1);
        ParentX                         = firstNetwork(counNodes,15):firstNetwork(counNodes,15)+firstNetwork(counNodes,18);
        ParentY                         = firstNetwork(counNodes,16):firstNetwork(counNodes,16)+firstNetwork(counNodes,19);
        ParentZ                         = firstNetwork(counNodes,17):firstNetwork(counNodes,17)+firstNetwork(counNodes,20);
        
        for counterChild =1:numChild
            %compare bounding boxes with children's
            currChildX                         = childLevel(counterChild,15):childLevel(counterChild,15)+childLevel(counterChild,18);
            currChildY                         = childLevel(counterChild,16):childLevel(counterChild,16)+childLevel(counterChild,19);
            currChildZ                         = childLevel(counterChild,17):childLevel(counterChild,17)+childLevel(counterChild,20);
            overlapInX                         = numel(intersect(ParentX,currChildX))/(numel(ParentX));
            overlapInY                         = numel(intersect(ParentY,currChildY))/(numel(ParentY));
            overlapInZ                         = numel(intersect(ParentZ,currChildZ))/(numel(ParentZ));
            
            % assign parenthood ONLY if there is > 30% overlap in X and Y and >60% in Z
            if (overlapInX>0.3)&&(overlapInY>0.3)&&(overlapInZ>0.6)
                % Assign parenthood ONLY if the child node has not got a parent yet
                if (childLevel(counterChild,7)==0)
                    %assing parent and child relationships
                    %---- parenthood is assigned to CHILD NODE
                    firstNetwork(childLevel(counterChild,6),7)    = firstNetwork(counNodes,6);
                    %---- childhood is assigned to  PARENT NODE
                    firstNetwork(counNodes,8)               = childLevel(counterChild,6);
                    %---- region is assigned to  -1 to describe that it was not with keyhole
                    firstNetwork(counNodes,12)               = -1;
                end
            end
        end
    end
end
