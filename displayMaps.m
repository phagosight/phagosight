function displayMaps(handles,dataR)
%function displayMaps (handles,dataR)
%
%--------------------------------------------------------------------------
%
% displayMaps  displays maps
%
%       INPUT
%         handles:      handles struct.
%         dataR:        path to reduced data.
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
figure


[rowsP,colsP,levsP]         = size(dataR);
if isstruct(handles)
    dataToMesh              = handles.distMaps.absDistMap;
    if isfield(handles,'ChannelDistribution')
        dataToSurf          = dataR(:,:,handles.ChannelDistribution(3));
    else
        dataToSurf          = dataR(:,:,levsP);       
    end
else
    dataToMesh              = handles;
    dataToSurf              = dataR(:,:,levsP);
end



meshdataR                   = zeros(rowsP,colsP);
[RR,CC]                     = meshgrid(1:colsP,1:rowsP);


dataToSurf                  = dataToSurf-min(dataToSurf(:));

stepMeshMap                 = 1;
stepMesh                    = 3;


hold off
handleS= surf(RR(1:stepMeshMap:end,1:stepMeshMap:end),CC(1:stepMeshMap:end,1:stepMeshMap:end),min(dataToMesh(:))+0.01+meshdataR(1:stepMeshMap:end,1:stepMeshMap:end), dataToSurf(1:stepMeshMap:end,1:stepMeshMap:end,1),'edgecolor','interp');
hold on

handleM = mesh(RR(1:stepMesh:end,1:stepMesh:end),CC(1:stepMesh:end,1:stepMesh:end),dataToMesh(1:stepMesh:end,1:stepMesh:end),'marker','none','LineStyle','-','facecolor','none','edgecolor',[0 0.8 0.8]);
axis tight
axis ij
set(handleM,'EdgeAlpha',0.4)
rotate3d on
view(45,40)
colormap(gray)
