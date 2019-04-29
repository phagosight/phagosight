function handles = joinMultipleTracks(handles,woundRegion,maxDistanceConsidered,minDistanceConsidered,minSeparationTwin)
%function handles = joinMultipleTracks(handles,woundRegion,maxDistanceConsidered,minDistanceConsidered,minSeparationTwin)
%
%
%--------------------------------------------------------------------------
% joinMultipleTracks analyses all tracks of handles, selects those that can be merged
%     under certain criteria and calls joinTwoTracks for every pair of candidates
%       Tasks:
%           1 measures distance from each track to upper tracks 
%           2 calculate distances if not provided by the user
%           3 discard pairs that do not satisfy criteria
%           4 call pairs into joinTwoTracks
%
%       INPUT
%         handles:                  handles struct
%         woundRegion:              wound region (to recalculate metrics)
%         maxDistanceConsidered:    above the maximum distance the pair is discarded
%         minDistanceConsidered:    below the minimum distance the pair is discarded
%         minSeparationTwin:        in case one track may be assigned as child to two
%                                   candidates, only assign if one is closer to the
%                                   child by this distance (default is to discard
%                                   both options and be done manually)
%
%       OUTPUT
%         handles:                  new handles struct with the tracks merged
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
numTracks = size(handles.finalNetwork,2);

if (~exist('woundRegion','var'))|(isempty(woundRegion))
    woundRegion = zeros(handles.rows,handles.cols);
    woundRegion (:,ceil(handles.cols/4):handles.cols) = 1;
end
%
%listTracks = 1:numTracks;

%%
recordTrackDistances =[];
for counterTrack = numTracks:-1:1  %39;
    currTrack                   = handles.nodeNetwork(handles.nodeNetwork(:,13)==counterTrack,1:5);
    timeIndex1                  = currTrack(:,5);
    currTrack_pos               = currTrack(end,1:5);
    max_diff                    = inf;
    max_diffT                   = inf;
    next_diff                   = inf;
    next_diffT                  = inf;
    candidateT                  = inf;
    %find closest initial track to final point of current track
    for counterTrack2           = counterTrack+1:numTracks
        %if listTracks(counterTrack2)>0
        nextTrack_pos       = handles.nodeNetwork(handles.finalNetwork(1,counterTrack2),1:5);
        timeIndex2          = handles.nodeNetwork(handles.nodeNetwork(:,13)==counterTrack2,5);
        diff_pos1            = sum(sqrt((currTrack_pos([1 2 3 ])-nextTrack_pos([1 2 3 ])).^2));
        diff_pos2            = (currTrack_pos([5 ])-nextTrack_pos([5 ]));
        diff_pos            = sqrt(diff_pos1.^2+(2*(abs(diff_pos2)+1)).^2);
        %discard joining if there is a total overlap of tracks
        if (numel(intersect(timeIndex2,timeIndex1))<min(numel(timeIndex1),numel(timeIndex2)))
            if diff_pos<max_diff
                if (abs(diff_pos2)-2)<abs(max_diffT)
                    
                    next_diff       = max_diff;
                    next_diffT      = max_diffT;
                    
                    max_diff        = diff_pos;
                    candidateT      = counterTrack2;
                    max_diffT        = diff_pos2;
                end
                
            end
        end
    end
    % The decision of joining two tracks is done on the basis of a distance between
    % end of one and beginning of the other. The distance is calculated on x,y,z,t so
    % the tracks can be more than 1 time frame appart from each other.
    recordTrackDistances=[recordTrackDistances;[round(max_diff) round(next_diff) round(max_diffT) counterTrack candidateT]];
    %disp([round(max_diff) round(next_diff) counterTrack candidateT])
end

%remove all cases in which the distance is INF
recordTrackDistances(recordTrackDistances(:,1)==inf,:)=[];

if (~exist('maxDistanceConsidered','var'))|(isempty(maxDistanceConsidered))
    maxDistanceConsidered = 2*std(recordTrackDistances(:,1))+mean(recordTrackDistances(:,1));
end

if (~exist('minDistanceConsidered','var'))|(isempty(minDistanceConsidered))
    minDistanceConsidered = mean(recordTrackDistances(:,1));
end


%remove all cases in which the distance is above a maximum, if maximum is not
%determined, it will be calculated as Mean+2*std
recordTrackDistances(recordTrackDistances(:,1)>maxDistanceConsidered,:)=[];


[q1,q2]     = hist(recordTrackDistances(:,5),unique(recordTrackDistances(:,5)));
q3          = q2(q1>1);
if ~exist('minSeparationTwin','var')
    % find duplicate cases (two tracks aim to merge with a single one)  and remove
    q4      = find(ismember(recordTrackDistances(:,5),q3));
    recordTrackDistances(q4,:)=[];
else
    % find duplicate cases and keep the closest if the difference is beyond a minimum
    for counterDuplic =1:numel(q3)
        indexDuplic =  find(recordTrackDistances(:,5)==q3(counterDuplic));
        if numel(indexDuplic)>2
            [differenceBetweenPairs1,indPairs]   = sort(recordTrackDistances(indexDuplic,1));

            recordTrackDistances(indexDuplic(indPairs(3:end)),:)=[];
            indexDuplic((indPairs(3:end)))=[];
            %indexDuplic = indexDuplic (indPairs(1:2));            
        end
        try
            differenceBetweenPairs   = diff(recordTrackDistances(indexDuplic,1));
        catch
            q=1;
        end
        if abs(differenceBetweenPairs) >minSeparationTwin
            % select the largest and remove
            if differenceBetweenPairs>0
                % remove the second case
                recordTrackDistances(indexDuplic(2),:)=[];
            else
                % remove the first case
                recordTrackDistances(indexDuplic(1),:)=[];
            end
        else
            % remove both cases and leave for manual analysis
            recordTrackDistances(indexDuplic,:)=[];
        end
    end
    
end

% discard those cases where the difference between closest and next is below a
% minimum distance
recordTrackDistances(:,6) = recordTrackDistances(:,2)-recordTrackDistances(:,1);
recordTrackDistances(recordTrackDistances(:,6)<minDistanceConsidered,:)=[];

numLinks = size(recordTrackDistances,1);
indexBeforeJoining = recordTrackDistances(:,5);
for counterLinks = 1:numLinks
    %disp(strcat('Processing link ',num2str(counterLinks),'/',num2str(numLinks)));
    candidate1              = recordTrackDistances(counterLinks,4);
    candidate2              = recordTrackDistances(counterLinks,5);
    numTracksPre            = size(handles.finalNetwork,2);
    handles                 = joinTwoTracks(handles,candidate1,candidate2,[]);
    numTracksPost           = size(handles.finalNetwork,2);
    if numTracksPre>numTracksPost
        indexBeforeJoining(indexBeforeJoining>candidate2) = indexBeforeJoining(indexBeforeJoining>candidate2)-1;
        recordTrackDistances(:,5) = indexBeforeJoining;
    end
 end
