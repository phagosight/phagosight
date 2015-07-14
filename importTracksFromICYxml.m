function handles = importTracksFromICYxml(filename)
%function handles = importTracksFromICYxml(filename)
%--------------------------------------------------------------------------
% This function translates the tracks from ICY xml format  to be
% readable phagosight. The tracks for ICY are organised as a list of tracks, each between
% labels like:
%       <track id="187512149"> </track>
% and inside there is one line per time point like this:
%       <detection
%       classname="plugins.nchenouard.particleTracking.sequenceGenerator.ProfileSpotTrack"
%       color="-20480" t="0" type="1" x="62.26605480706269" y="12.835458913910855"
%       z="0"/>
% therefore the translation is rather more complex than from phagosight to ICY. In
% addition some parameters like volume or distance to closest neutrophil are set to
% zero.
% It is necessary to loop several times to generate the fields nodeNetwork and
% finalNetwork.
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
% This m-file is part of the PhagoSight package used to analyse
% fluorescent phagocytes as observed through confocal or
% multiphoton microscopes.  For a comprehensive 
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
%filename = 'HK_mat_Ha/tracks2.xml';
if ~exist('filename','var')
    filename = 'tracks.xml';
end

% First, open the file
f = fopen( filename, 'r' ); 

%%
% read the first IDs
line = fgets(f); % fprintf(f, '%s\n','<?xml version="1.0" encoding="UTF-8" standalone="no"?>');
line = fgets(f);% fprintf(f, '%s\n','<root>');
line = fgets(f);% fprintf(f, '%s\n','<trackfile version="1"/>');
line = fgets(f);% fprintf(f, '%s\n','<trackgroup description="refTracks0">');

% the string for each position will be based on a basic string split in parts
% str1 = '<detection classname="plugins.nchenouard.particleTracking.sequenceGenerator.ProfileSpotTrack" color="-20480" t="';
% str2 = '" type="1" x="';
% str3 = '" y="';
% str4 = '" z="';
% str5 = '"/>';
%%
% The tracks will be read one by one, and then each track will be appended to a
% struct with one time point per cell and each object will be one line, at the end,
% nodeNetwork can be formed by appending the times of the cells into a single matrix.

tracksPerTrack  = {};
tracksTime      = {};


trackID         = 0 ;
%%
% read now until </trackgroup> is found

while isempty(strfind(line,'</trackgroup>'))
    line = fgets(f);
    startTrack          = strfind(line,'<track id="');
    endTrack            = strfind(line,'</track');
    if ~isempty(startTrack)
        % this is a track that starts
        delimiterTrack  = strfind(line,'">');
        %trackID         = str2num (line(12:delimiterTrack-1));
        trackID         = trackID + 1;
        % initialise the cell for the current track
        tracksPerTrack{trackID}=[];
    end
    if (isempty(startTrack))&(isempty(endTrack))
        % this is a node of the track, analyse
        timePositionS   = strfind(line,'t="');
        timePositionE   = strfind(line,'" type="');
        xPosition       = strfind(line,'x="');
        yPosition       = strfind(line,'y="');
        zPosition       = strfind(line,'z="');     
        
        time            = 1+str2num (line(timePositionS+3:timePositionE-1));
        x               = str2num (line(xPosition+3:yPosition-3));
        y               = str2num (line(yPosition+3:zPosition-3));
        z               = str2num (line(zPosition+3:end-4));
        currObject      = [x y z 0 time zeros(1,7) trackID zeros(1,19)];
           
        tracksPerTrack{trackID}=[tracksPerTrack{trackID};currObject];    
        if size(tracksTime,1)<time
            tracksTime{time,1}=[];
        end

        tracksTime{time,1}=[tracksTime{time,1};currObject];    

    end
        % this is not a track, it may be a line or the end of the track

end
%% assign unique IDs into tracksTime
uniqueID =0;
for counterTime =1:size(tracksTime,1)
    uniqueID_end                    = uniqueID+size(tracksTime{counterTime},1);
    tracksTime{counterTime}(:,6)    = (uniqueID+1:uniqueID+size(tracksTime{counterTime},1));
    uniqueID                        = uniqueID_end;
end
%% assign parenthood and create finalNetwork
finalNetwork =[];
for counterTrack =1:size(tracksPerTrack,2)
    for counterTimeInTrack =1:size(tracksPerTrack{counterTrack},1)-1
        time_1  = tracksPerTrack{counterTrack}(counterTimeInTrack,5);
        time_2  = tracksPerTrack{counterTrack}(counterTimeInTrack+1,5);
        %parent
        indexParent             = find(tracksTime{time_1}(:,13)==counterTrack);
 
        indexChild              = find(tracksTime{time_2}(:,13)==counterTrack);

        currentParent           = tracksTime{time_1}(indexParent,:);  
        currentChild            = tracksTime{time_2}(indexChild,:);  

        tracksTime{time_1}(indexParent,8)    =currentChild(6);

        tracksTime{time_2}(indexChild,7)     =currentParent(6);

        finalNetwork (counterTimeInTrack:counterTimeInTrack+1,counterTrack) = [currentParent(6);currentChild(6)];

    end
end
%% Re arrage into nodeNetwork 
nodeNetwork =[];
for counterTime =1:size(tracksTime,1)
    nodeNetwork =[nodeNetwork;tracksTime{counterTime}];
end
% assign unique ID
nodeNetwork(:,6)        = (1:size(nodeNetwork,1));


handles.numFrames       = nodeNetwork(end,5);
handles.rows            = ceil(max(nodeNetwork(:,1)));
handles.cols            = ceil(max(nodeNetwork(:,2)));
handles.levs            = ceil(max(nodeNetwork(:,3)));
handles.nodeNetwork     = nodeNetwork;
handles.finalNetwork    = finalNetwork;



