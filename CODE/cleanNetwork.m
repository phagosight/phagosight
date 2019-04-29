function reducedNetwork = cleanNetwork(finalNetwork,firstNetwork)
%function reducedNetwork = cleanNetwork(finalNetwork,firstNetwork)
%
%--------------------------------------------------------------------------
% reduceNetwork  cleans the network from 2-node branches
%
%       INPUT
%           finalNetwork:   finalNewtwork
%           firstNetwork:   firstNewtwork
%       OUTPUT
%           reducedNetwork:	reduced network after cleaning
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
%     Foobar is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
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

    trackLength                     = sum(finalNetwork>0);
    reducedNetwork                 = finalNetwork(:,trackLength>0);

    %----- check the FIRST hop, since this is the most uncertain one.
    trackLength                     = sum(reducedNetwork>0);
    for counterTrack=1:size(reducedNetwork,2)
        if trackLength(counterTrack)>2
            %----- test the FIRST node of the Track
            childLevel                  = firstNetwork(reducedNetwork(1,counterTrack),:);
            %----- against the SECOND and THIRD nodes of the same branch
            keyHoleDef                  = [firstNetwork(reducedNetwork(2,counterTrack),1:3); firstNetwork(reducedNetwork(3,counterTrack),1:3) ];

            inTheTail               = testKeyHole(keyHoleDef,childLevel(:,1:3),firstNetwork(reducedNetwork(2,counterTrack) ,4));
            
            %------ if there is no node in the tail delete the first node of the track
            if (inTheTail==0)
                reducedNetwork(1:trackLength(counterTrack)-1,counterTrack)=reducedNetwork(2:trackLength(counterTrack),counterTrack);
                reducedNetwork(trackLength(counterTrack),counterTrack)=0;
            end
        end
    end
end
