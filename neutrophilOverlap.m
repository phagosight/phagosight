function [finalLabel] = neutrophilOverlap(firstNetwork,linkedNetwork,handles)
%function [finalLabel] = neutrophilOverlap(firstNetwork,linkedNetwork,handles)
% 
%
%--------------------------------------------------------------------------
% neutrophilOverlap  analyses if neutrophils overlap between time frames
%
%       INPUT
%         firstNetwork:         first network matrix
%         linkedNetwork:        linked network matrix
%         handles:              handles struct
%
%       OUTPUT
%         finalLabel:           labels for tracks
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


%%
beginTrackNode                                              = linkedNetwork(1,:);
endTrackNode                                                = max(linkedNetwork);

beginTrackFrame                                             = firstNetwork(beginTrackNode,5);
endTrackFrame                                               = firstNetwork(endTrackNode,5);


possibleMergerT                                             = unique(endTrackFrame);
%%
%discard tracks that finish at numTracks
possibleMergerT(possibleMergerT==handles.numFrames)         = [];

mergeTrackRes                                               = [];
splitTrackRes                                               = [];
%------------ first test the mergers
if ~isempty(possibleMergerT)
    for counterMerge = 1:size(possibleMergerT,1)
        %for all frames where a track finishes, check if there are at least two tracks that finish
        if (sum(endTrackFrame==possibleMergerT(counterMerge))==2)
            %now check if there is a track that start immediately after that
            if sum(beginTrackFrame==(1+possibleMergerT(counterMerge)))>0
                %Determine the objects of interest
                mergingTrack                                    = find(endTrackFrame==possibleMergerT(counterMerge));
                mergedTrack                                     = find(beginTrackFrame==(1+possibleMergerT(counterMerge)));
                Object1                                         = endTrackNode(mergingTrack(1));
                Object2                                         = endTrackNode(mergingTrack(2));
                for countMergeOptions=1:numel(mergedTrack)
                    NewObject                                       = beginTrackNode(mergedTrack(countMergeOptions));
                    %check that the centroids/bodies of the merging objects overlap with the merged one
                    xBox                                            = [firstNetwork([Object1 Object2 NewObject],15) sum(firstNetwork([Object1 Object2 NewObject],[15 18]),2)];
                    yBox                                            = [firstNetwork([Object1 Object2 NewObject],16) sum(firstNetwork([Object1 Object2 NewObject],[16 19]),2)];
                    zBox                                            = [firstNetwork([Object1 Object2 NewObject],17) sum(firstNetwork([Object1 Object2 NewObject],[17 20]),2)];

                    condX1                                          = sum(ismember(xBox(1,1):xBox(1,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(1,1)-xBox(1,2)));
                    condY1                                          = sum(ismember(yBox(1,1):yBox(1,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(1,1)-yBox(1,2)));
                    condZ1                                          = sum(ismember(zBox(1,1):zBox(1,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(1,1)-zBox(1,2)));

                    condX2                                          = sum(ismember(xBox(2,1):xBox(2,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(2,1)-xBox(2,2)));
                    condY2                                          = sum(ismember(yBox(2,1):yBox(2,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(2,1)-yBox(2,2)));
                    condZ2                                          = sum(ismember(zBox(2,1):zBox(2,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(2,1)-zBox(2,2)));

                    overLapFut(countMergeOptions)                   = condX1*condY1*condZ1*condX2*condY2*condZ2; %#ok<AGROW>
                    

                end
                if  any(overLapFut>0.5)
                    [q1,q2] = max(overLapFut);
                    
                    mergeTrackRes                               = [mergeTrackRes;mergingTrack' mergedTrack(q2)]; %#ok<AGROW>
                    
                end
                clear overLapFut*
            end
        end
    end



    %------------ second test the split
    for counterSplit = 1:size(possibleMergerT,1)
        %for all frames where a track finishes, check if there are at least two tracks that begin after
        if (sum(endTrackFrame==possibleMergerT(counterSplit))>0)
            %now check if there is a track that start immediately after that
            if sum(beginTrackFrame==(1+possibleMergerT(counterSplit)))>1
                %Determine the objects of interest
                splittingTrack                                  = find(endTrackFrame==possibleMergerT(counterSplit));
                splittedTrack                                   = find(beginTrackFrame==(1+possibleMergerT(counterSplit)));
                OldObject                                       = endTrackNode(splittingTrack(1));
                numSplitOptions                                 = numel(splittedTrack);

                if (numSplitOptions ==2)

                    Object1                                         = beginTrackNode(splittedTrack(1));
                    Object2                                         = beginTrackNode(splittedTrack(2));

                    %check that the centroids/bodies of the merging objects overlap with the merged one
                    xBox                                            = [firstNetwork([Object1 Object2 OldObject],15) sum(firstNetwork([Object1 Object2 OldObject],[15 18]),2)];
                    yBox                                            = [firstNetwork([Object1 Object2 OldObject],16) sum(firstNetwork([Object1 Object2 OldObject],[16 19]),2)];
                    zBox                                            = [firstNetwork([Object1 Object2 OldObject],17) sum(firstNetwork([Object1 Object2 OldObject],[17 20]),2)];

                    condX1                                          = sum(ismember(xBox(1,1):xBox(1,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(1,1)-xBox(1,2)));
                    condY1                                          = sum(ismember(yBox(1,1):yBox(1,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(1,1)-yBox(1,2)));
                    condZ1                                          = sum(ismember(zBox(1,1):zBox(1,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(1,1)-zBox(1,2)));

                    condX2                                          = sum(ismember(xBox(2,1):xBox(2,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(2,1)-xBox(2,2)));
                    condY2                                          = sum(ismember(yBox(2,1):yBox(2,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(2,1)-yBox(2,2)));
                    condZ2                                          = sum(ismember(zBox(2,1):zBox(2,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(2,1)-zBox(2,2)));

                    overLapFut                                      = condX1*condY1*condZ1*condX2*condY2*condZ2;
                    if  overLapFut>0.5
                        
                        splitTrackRes                               = [splitTrackRes;splittingTrack(1) splittedTrack']; %#ok<AGROW>
                        
                    end
                    clear overLapFut*

                else
                    for count1SplitOptions = 1: numSplitOptions -1
                        for count2SplitOptions = 1+count1SplitOptions : numSplitOptions

                            Object1                                         = beginTrackNode(splittedTrack(count1SplitOptions));
                            Object2                                         = beginTrackNode(splittedTrack(count2SplitOptions));

                            %check that the centroids/bodies of the merging objects overlap with the merged one
                            xBox                                            = [firstNetwork([Object1 Object2 OldObject],15) sum(firstNetwork([Object1 Object2 OldObject],[15 18]),2)];
                            yBox                                            = [firstNetwork([Object1 Object2 OldObject],16) sum(firstNetwork([Object1 Object2 OldObject],[16 19]),2)];
                            zBox                                            = [firstNetwork([Object1 Object2 OldObject],17) sum(firstNetwork([Object1 Object2 OldObject],[17 20]),2)];

                            condX1                                          = sum(ismember(xBox(1,1):xBox(1,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(1,1)-xBox(1,2)));
                            condY1                                          = sum(ismember(yBox(1,1):yBox(1,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(1,1)-yBox(1,2)));
                            condZ1                                          = sum(ismember(zBox(1,1):zBox(1,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(1,1)-zBox(1,2)));

                            condX2                                          = sum(ismember(xBox(2,1):xBox(2,2),xBox(3,1):xBox(3,2)))/(1+abs(xBox(2,1)-xBox(2,2)));
                            condY2                                          = sum(ismember(yBox(2,1):yBox(2,2),yBox(3,1):yBox(3,2)))/(1+abs(yBox(2,1)-yBox(2,2)));
                            condZ2                                          = sum(ismember(zBox(2,1):zBox(2,2),zBox(3,1):zBox(3,2)))/(1+abs(zBox(2,1)-zBox(2,2)));

                            overLapFut(count1SplitOptions,count2SplitOptions)= condX1*condY1*condZ1*condX2*condY2*condZ2;
                            overLapFut1(count1SplitOptions,count2SplitOptions)=count1SplitOptions; %#ok<AGROW>
                            overLapFut2(count1SplitOptions,count2SplitOptions)=count2SplitOptions; %#ok<AGROW>
                        end
                    end
                    if  any(overLapFut(:)>0.5)
                        maxOverL                                        = max(overLapFut(:));
                        maxOverCoord                                    = find(overLapFut==maxOverL,1);

                        %try
                        splitTrackRes                               = [splitTrackRes;splittingTrack(1) splittedTrack(overLapFut1(maxOverCoord)) splittedTrack(overLapFut2(maxOverCoord))]; %#ok<AGROW>
                        %catch
                        %    qqq=1;
                        %end
                    end
                    clear overLapFut*

                end
            end
        end
    end

    % mergeTrackRes     - contains the tracks that merge into a single track    [Track1 Track2 NewTrack]
    % splitTrackRes     - contains the tracks that split into two track         [OldTrack Track1 Track2]

    % tt is a matrix that keeps track of merge and splits, row 4 contains all tracks.
    numTracks                                                   = size(linkedNetwork,2);
    tt                                                          = zeros(6,numTracks);
    tt(4,:)                                                     = 1:numTracks;

    % loop to obtain average, initial and final volume of tracks
    for k=1:numTracks
        lenTrack{k}                                             = sum(linkedNetwork(:,k)>0); %#ok<AGROW>
        volTrack{k}                                             = firstNetwork(linkedNetwork(1:lenTrack{k},k),10);%#ok<AGROW>
        avVolTrack(1,k)                                         = round(mean(volTrack{k}));%#ok<AGROW>
        initVolTrack(1,k)                                       = round((volTrack{k}(1)));%#ok<NASGU,AGROW>
        finVolTrack(1,k)                                        = round((volTrack{k}(end)));%#ok<NASGU,AGROW>
    end

    % Below line 4, for each track that is merged the new track (ONE number at row 5) number will be placed
    for k=1:size(mergeTrackRes,1)
        tt(5,mergeTrackRes(k,1))                                = mergeTrackRes(k,3);
        tt(5,mergeTrackRes(k,2))                                = mergeTrackRes(k,3);
    end


    % when lines are split, the new numbers (TWO one at row 5 and one at 6) are is placed below
    for k=1:size(splitTrackRes,1)
        tt(5,splitTrackRes(k,1))                                = splitTrackRes(k,2);
        tt(6,splitTrackRes(k,1))                                = splitTrackRes(k,3);
    end


    % This loop will assign to those tracks that start from a split, its original parent
    for k=1:numTracks

        if tt(5,k)~=0
            %either a parent of a merged node or the parent of 2 split nodes
            if tt(6,k)==0
                %parent of merger
                if tt(1,tt(5,k))==0
                    %first location has not been taken
                    if tt(1,k)~=0           %assign the parenthood (FIRST LINE) as the parent the track has already been assigned
                        tt(1,tt(5,k))                           = (tt(1,k) );
                    else                    %assign the parenthood (FIRST LINE) as its own number track
                        tt(1,tt(5,k))                           = k;
                    end
                    %tt(1,tt(5,k))=max(k,tt(1,k) );
                else
                    %first location taken, pass to second
                    if tt(1,k)~=0           %assign the parenthood (SECOND LINE) as the parent the track has already been assigned
                        tt(2,tt(5,k))                           = (tt(1,k) );
                    else                    %assign the parenthood (SECOND LINE) as its own number track
                        tt(2,tt(5,k))                           = k;
                    end

                end
            else
                %parent of split
                if (tt(1,k)~=0)&(tt(2,k)~=0)               
                    %previous parent assigned in row 1 and 2
                    %in order to decide which parent is assigned to each child a volume comparison has to be
                    %performed, compare areas 1-3 2-4 and 1-4 2-3 to decide how are the children assigned
                        firstComp                                   = sum( [(avVolTrack(tt(1,k))-avVolTrack(tt(5,k))) (avVolTrack(tt(2,k))-avVolTrack(tt(6,k)))].^2 );
                        secondComp                                  = sum( [(avVolTrack(tt(2,k))-avVolTrack(tt(5,k))) (avVolTrack(tt(1,k))-avVolTrack(tt(6,k)))].^2 );
                    
                    %once parenthood is assigned, the volumes should also be updated as neutrophils may be
                    %losing volume (going up or down the field of view or losing intensity and therefore
                    %volume)
                    if  firstComp<secondComp
                        tt(1,tt(5,k))                           = (tt(1,k) );
                        tt(1,tt(6,k))                           = (tt(2,k) );
                        % average volume values
                        avVolTrack(tt(1,k))                     = 0.5*(avVolTrack(tt(1,k))+avVolTrack(tt(5,k))); %#ok<AGROW>
                        avVolTrack(tt(2,k))                     = 0.5*(avVolTrack(tt(2,k))+avVolTrack(tt(6,k))); %#ok<AGROW>
                    else
                        tt(1,tt(5,k))                           = (tt(2,k) );
                        tt(1,tt(6,k))                           = (tt(1,k) );
                        % average volume values
                        avVolTrack(tt(2,k))                     = 0.5*(avVolTrack(tt(2,k))+avVolTrack(tt(5,k))); %#ok<AGROW>
                        avVolTrack(tt(1,k))                     = 0.5*(avVolTrack(tt(1,k))+avVolTrack(tt(6,k))); %#ok<AGROW>
                    end
                elseif (tt(1,k)~=0)&(tt(2,k)==0)               
                    %previous parent assigned in row 1 and 2 hs no parent, it is a cell that comes from a subdivision,
                    %assign common parent
                    tt(1,tt(5,k))                               = tt(1,k);
                    tt(1,tt(6,k))                               = tt(1,k);
                else
                    %No previous parent
                    tt(1,tt(5,k))                               = k;
                    tt(1,tt(6,k))                               = k;

                end

            end

        end

    end
    tt(tt==0)                                                   = inf;
    finalLabel                                                  = min(tt([1 4],:));
else
    numTracks                                                   = size(linkedNetwork,2);
    finalLabel                                                  = 1:numTracks;
end


