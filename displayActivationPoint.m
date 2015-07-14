function displayActivationPoint (handles,dataR,woundRegion)
%function displayActivationPoint (handles,dataR,woundRegion)
%
%--------------------------------------------------------------------------
%
% displayActivationPoint  displays the point where the neutrophil is activated 
%     as calculated previously according to its change of slope with time
%
%       INPUT
%         handles:      handles structure.
%         dataR:        either a volume or a slice of the DIC data to be displayed with the tracks
%         woundRegion:  the region that has been manually drawn to define the region of the wound  
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

if (~isfield(handles,'distMaps'))||(~isfield(handles.distMaps,'activationPoint3'))
    if ~exist('woundRegion','var')
        disp('This function requires the parameters calculated with effectiveDistance and effectiveTracks2.');
        disp('You can obtain them before calling this function, or pass the variable woundRegion as an input argument');
        return;
    else
        handles                                     = effectiveDistance(handles,woundRegion);
        handles                                     = effectiveTracks(handles,woundRegion);
    end
end

hold off

%plotTracks(handles.nodeNetwork(:,[1 2 5]), handles.finalNetwork(:,handles.distanceNetwork.numHops>0 ),2,1,1);axis([1 handles.rows 1 handles.cols 1 handles.numFrames])
plotTracks(handles,2);axis([1 handles.cols 1 handles.rows 1 handles.numFrames])
hold on

indexP                                          = find(handles.distMaps.activationTime3~=0);
plot3(handles.distMaps.activationPoint3(1,indexP),handles.distMaps.activationPoint3(2,indexP),handles.distMaps.activationTime3(indexP),'ko','markersize',9)
indexP                                          = find(handles.distMaps.enterWoundTime3~=0);

plot3(handles.distMaps.enterWoundPoint3(1,indexP),handles.distMaps.enterWoundPoint3(2,indexP),handles.distMaps.enterWoundTime3(indexP),'r*','markersize',9)


view(0,0)

%display the DIC channel with the tracks. Select the channel of the DIC
if exist('dataR','var')
    if size(dataR,3)==1
        dataR_toPlot    = dataR;
    else
        if handles.ChannelDistribution(3)==0
            disp('The handles.ChannelDistribution does not have a DIC channel indicated.');
            return;
        else
            dataR_toPlot    = dataR(:,:,handles.ChannelDistribution(3));
        end
    end
    [rows,cols,levs]    = size(dataR);
    [YY,XX]             = meshgrid(1:cols,1:rows);
    ZZ                  = zeros(rows,cols);
    h1                  = surf(YY,XX,ZZ,dataR_toPlot,'edgecolor','none');
    colorbar off
    allObjects          = findall(gcf);
    hh=findobj(allObjects,'tag','colorbarLabel');
    set(hh(1),'visible','off')
    colormap(gray)
end
    
    
