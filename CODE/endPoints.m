function [EndPoints1,EndPoints2,numEndPoints]=endPoints(dataIn)
%function [EndPoints1,EndPoints2,numEndPoints]=endPoints(dataIn)
%
%--------------------------------------------------------------------------
% endPoints  locates the end points of lines of skeletonised structures
%
%       INPUT
%         dataIn:           data with lines on it.
%
%       OUTPUT
%         EndPoints1:       1 where the lines finish, 0 elsewhere
%         EndPoints2:       Endpoints1 dilated for better visual display
%         numEndPoints:     number of end points detected
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

EndKernel1                   = [ 0 0 0; -1 1 -1; -1 -1 -1];    %End points
EndKernelEx                  = [ 1 0 1; -1 1 -1; -1 -1 -1];    %End points

dataIn=(padData(dataIn,1,[],0));

EndPoints1                      =zeros(size(dataIn));
EndPointsEx                     = zeros(size(dataIn));
%%
for k=0:3
    EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel1,k)));
end

EndPoints1                   = EndPoints1(2:end-1,2:end-1);
%%


for k=0:3
    EndPointsEx               = (EndPointsEx|bwhitmiss(dataIn,rot90(EndKernelEx,k)));
end

EndPointsEx                   = EndPointsEx(2:end-1,2:end-1);

%%
EndPoints1 = EndPoints1 &(~EndPointsEx);

%%
if nargout>=2
    EndPoints2               = EndPoints1;
    EndPoints2(1:end-1,:)    = (EndPoints1(1:end-1,:)|EndPoints1(2:end,:));
    EndPoints2(2:end,:)      = (EndPoints1(1:end-1,:)|EndPoints2(2:end,:));
    EndPoints2(:,1:end-1)    = (EndPoints2(:,1:end-1)|EndPoints1(:,2:end));
    EndPoints2(:,2:end)      = (EndPoints1(:,1:end-1)|EndPoints2(:,2:end));
    
    numEndPoints             = sum(find(EndPoints1)>0);
    
end
