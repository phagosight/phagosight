function distanceNetwork = getDistanceNetForTracks(finalNetwork,parentNetwork)
%%function distanceNetwork = = getDistanceNetForTracks(finalNetwork,parentNetwork)
%
%--------------------------------------------------------------------------
% distanceNetwork  computes the distance network for the nodes contained in
%     the final and parent networks
%
%       INPUT
%         finalNetwork:     final network
%         parentNetwork:    parent network
%
%       OUTPUT
%         distanceNetowrk:  distance network
%
%     This function used to be in linkEndStartTracks.m under the name
%     getDistanceNet, as that name was in conflict with the getDistanceNet
%     extracted from trackNeutrophil, the name was changed to
%     getDistanceNetForTracks
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


    numNodesBranches                = sum(finalNetwork>0);
    for k=1:size(finalNetwork,2)
        indexN=finalNetwork(:,k);indexN(indexN==0)=[];
        %numNodes=size(indexN,1);
        distPerBranch(1:numNodesBranches(k)-1,k)=sqrt(sum((diff(parentNetwork(indexN,1:3))).^2,2)); %#ok<AGROW>
    end

    distanceNetwork.perHop          = distPerBranch;
    distanceNetwork.numHops         = numNodesBranches;
    distanceNetwork.totPerTrack     = sum(distPerBranch);
    distanceNetwork.maxPerTrack     = max(distPerBranch);
    distanceNetwork.avPerTrack      = distanceNetwork.totPerTrack./(numNodesBranches-1);
end
