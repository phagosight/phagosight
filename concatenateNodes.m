function neighbourNetwork=concatenateNodes(parentNetwork)
%function neighbourNetwork=concatenateNodes(parentNetwork)
%
%--------------------------------------------------------------------------    
% concatenateNodes  creates a matrix containing the neighbourhood
%     relationship between nodes
%
%       INPUT
%           parentNetwork:      parent Network
%       OUTPUT
%           neighbourNetwork:   matrix with information about node's
%                               neighbourhood
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

    
    
    %------ can we get distances here?
    [numNodesInNet] = size(parentNetwork,1);
    counNewBranches = 1;
    %----- RE-assign an ID for the nodes of the Network (some may be removed later on, so keep this to identify parent/child
    parentNetwork(:,6)=(1:numNodesInNet)';

    neighbourNetwork=[];
    %------ concatenate parent/child nodes
    for counNodes=1:numNodesInNet
        %------ will loop over all the nodes, check if they have been assigned to a branch by
        %------ looking at their ID which will be erased in case it is and do a second loop over their children
        %------ the structure of parentNetwork (x,y,z,d, timeSlice-5, ID-6, Parent-7, child-8)
        presentNode=counNodes;

        if parentNetwork(presentNode,6)~=0
            %------ the node has NOT been assigned to a branch, it should be assigned to a new Branch
            if parentNetwork(presentNode,8)~=0
                %------ the node has a child  then proceed
                childNode=parentNetwork(presentNode,8);
                if parentNetwork(childNode,6)~=0
                    %------ child node has not been assigned yet, then assign
                    counBranchDepth=1;
                    neighbourNetwork(counBranchDepth,counNewBranches)=parentNetwork(presentNode,6); %#ok<AGROW>
                    counBranchDepth=        counBranchDepth+1;
                    neighbourNetwork(counBranchDepth,counNewBranches)=parentNetwork(childNode,6); %#ok<AGROW>
                    %------- erase IDs of parent and child
                    parentNetwork(presentNode,6)=0;
                    parentNetwork(childNode,6)=0;

                    while parentNetwork(childNode,8)~=0
                        %------ there is a child, continue process downwards
                        presentNode=childNode;
                        childNode=parentNetwork(presentNode,8);
                        if parentNetwork(childNode,6)~=0
                            %------ child node has not been assigned yet, then assign
                            counBranchDepth=        counBranchDepth+1;
                            neighbourNetwork(counBranchDepth,counNewBranches)=parentNetwork(presentNode,8); %#ok<AGROW>
                            %------- erase IDs of parent and child
                            parentNetwork(presentNode,6)=0;
                            parentNetwork(childNode,6)=0;
                        end
                    end
                    counNewBranches=counNewBranches+1;
                end
            end
        end
    end
    %-------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------
end
