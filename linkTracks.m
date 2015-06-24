function linkedNetwork = linkTracks(reducedNetwork,firstNetwork)
%linkedNetwork = linkTracks(reducedNetwork,firstNetwork)
%
%--------------------------------------------------------------------------
% linkTracks uses the reducedNetwork and firstNetwork to link tracks
%     based in the keyhole
%
%       INPUT
%         reducedNetwork:           reduced network matrix (included in handles)
%         firstNetwork:             first network matrix
%
%       OUTPUT
%         linkedNetwork:            updated network matrix with linked tracks
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
    %-----
    [longestTrack,numTracks]        = size(reducedNetwork);
    trackLength                     = sum(reducedNetwork>0);
    linkedNetwork=reducedNetwork;
    %----- backtrack the tacks to see if they can be linked together
    tracksToLink=[];
    for counterStart=1:numTracks
        if linkedNetwork(1,counterStart)>0
            vectorOfTarget=[1:counterStart-1 counterStart+1:numTracks ];
            for counterTarget=vectorOfTarget
                %        for counterTarget=counterStart+1:numTracks
                if linkedNetwork(1,counterTarget)>0
                    %----- test the LAST node of the Track
                    childLevel                  = firstNetwork(linkedNetwork(trackLength(counterStart)  ,counterStart),:);
                    childFrame                  = firstNetwork(linkedNetwork(trackLength(counterStart)  ,counterStart),5);
                    %----- against the FIRST and SECOND nodes of the other tracks
                    %                keyHoleDef                  = [firstNetwork(linkedNetwork(1,counterTarget),1:3); firstNetwork(linkedNetwork(2,counterTarget),1:3) ];
                    %----- against the WHOLE BRANCH of the other tracks
                    keyHoleDef                  = (firstNetwork(linkedNetwork(1:trackLength(counterTarget),counterTarget),1:3) );
                    keyHoleFrame                = (firstNetwork(linkedNetwork(1,counterTarget),5));

                    if childFrame==(keyHoleFrame-1)

                        inTheTail               = testKeyHole(keyHoleDef,childLevel(:,1:3),firstNetwork(linkedNetwork(2,counterTarget) ,4));
                        
                        %------ if there is no node in the tail delete the first node of the track
                        if (inTheTail==0)
                            %------ no connection between branches has been detected
                        else
                            %------ the branches have a possible connection, keep record and act after all
                            %------ branches have been analysed
                            tracksToLink=[tracksToLink;[counterStart counterTarget]]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end



   
    %%
    %------ check that there is no branch that can link to two different branches either above or below
    numBranchesToLink                           = size(tracksToLink,1);
    for counterStart = 1:numBranchesToLink-1
        tempStart                               = tracksToLink(counterStart,1);
        tempTarget                              = tracksToLink(counterStart,2);
        detectChange=0;
        for counterTarget = counterStart+1:numBranchesToLink
            if (tracksToLink(counterTarget,1)==tempStart)||(tracksToLink(counterTarget,2)==tempTarget)
                tracksToLink(counterTarget,:)   = 0;
                detectChange=1;

            end
        end
        if detectChange==1
            tracksToLink(counterStart,:)    = 0;
        end
    end
    %------- eliminate the branches
    if ~isempty(tracksToLink)
        indexTracks=find(tracksToLink(:,1));
        if sum(indexTracks>0)
            tracksToLink=tracksToLink(indexTracks,:);
            numBranchesToLink                           = size(tracksToLink,1);

            for counterTrack=numBranchesToLink:-1:1
                secondBranch                            = linkedNetwork(:,tracksToLink(counterTrack,2));
                sizeSecondBranch                        = sum(secondBranch>0);
                firstBranch                             = linkedNetwork(:,tracksToLink(counterTrack,1));
                sizeFirstBranch                         = sum(firstBranch>0);
                tempSize                                = sizeFirstBranch+sizeSecondBranch;
                linkedNetwork(tempSize,tracksToLink(counterTrack,1))=0;
                linkedNetwork(sizeFirstBranch+1:tempSize,tracksToLink(counterTrack,1))=secondBranch(1:sizeSecondBranch);
                linkedNetwork(:,tracksToLink(counterTrack,2))=0;
                    trackLength                     = sum(linkedNetwork>0);
            end
        end
    end
    trackLength                                 = sum(linkedNetwork>0);
    linkedNetwork                               = linkedNetwork(:,trackLength>0);
end
