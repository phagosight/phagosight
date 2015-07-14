function [BranchPoints1,BranchPoints2,numBranchPoints] = branchPoints(dataIn)
%function [BranchPoints1,BranchPoints2,numBranchPoints] = branchPoints(dataIn)
%
%--------------------------------------------------------------------------
% branchPoints  locates the branching points of lines of skeletonised 
%     structures
%
%       INPUT
%           dataIn:             data with lines on it.
%
%       OUTPUT  
%           BranchPoints1:      1 where the lines divide, 0 elsewhere.
%           BranchPoints2:      endpoints1 dilated for better visual display.
%           numBranchPoints:    number of branch points detected.
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

BranchKernel1                   = [0 1 0; 1 1 0; 0 0 1];    %Corner
BranchKernel2                   = [0 1 0; 0 1 0; 1 0 1];    %Y
BranchKernel3                   = [1 0 0; 0 1 0; 1 0 1];    %diagonal
%BranchKernel4                  = [1 1 1; 0 1 0; 0 1 0];    %T

dataIn=(padData(dataIn,1,[],0));

BranchPoints1=zeros(size(dataIn));

for k=0:3
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel1,k)));
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel2,k)));
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel3,k)));
end

BranchPoints1                   = BranchPoints1(2:end-1,2:end-1);



if nargout>=2
    BranchPoints2               = BranchPoints1;
    BranchPoints2(1:end-1,:)    = (BranchPoints1(1:end-1,:)|BranchPoints1(2:end,:));
    BranchPoints2(2:end,:)      = (BranchPoints1(1:end-1,:)|BranchPoints2(2:end,:));
    BranchPoints2(:,1:end-1)    = (BranchPoints2(:,1:end-1)|BranchPoints1(:,2:end));
    BranchPoints2(:,2:end)      = (BranchPoints1(:,1:end-1)|BranchPoints2(:,2:end));
    
    numBranchPoints             = sum(find(BranchPoints1)>0);
    
end
