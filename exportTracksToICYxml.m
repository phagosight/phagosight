function exportTracksToICYxml(handles,filename)

% translate the tracks from phagosight format to an xml file to be readable by ICY

% The tracks for ICY are organised as a list of tracks, each between labels like:
%       <track id="187512149">
%       </track>
% and inside there is one line per time point like this:
%       <detection classname="plugins.nchenouard.particleTracking.sequenceGenerator.ProfileSpotTrack" color="-20480" t="0" type="1" x="62.26605480706269" y="12.835458913910855" z="0"/>
% therefore it is necessary to loop over handles.finalNetwork to create each track,
% and then loop over the points of the track.

if ~exist('filename','var')
    filename = 'tracks.xml';
end
% First, open the file
f = fopen( filename, 'w' ); 
% print the first IDs
fprintf(f, '%s\n','<?xml version="1.0" encoding="UTF-8" standalone="no"?>');
fprintf(f, '%s\n','<root>');
fprintf(f, '%s\n','<trackfile version="1"/>');
fprintf(f, '%s\n','<trackgroup description="refTracks0">');

% the string for each position will be based on a basic string split in parts
str1 = '<detection classname="plugins.nchenouard.particleTracking.sequenceGenerator.ProfileSpotTrack" color="-20480" t="';
str2 = '" type="1" x="';
str3 = '" y="';
str4 = '" z="';
str5 = '"/>';


for counterTrack=1:size(handles.finalNetwork,2)
    % start the track, add labels
    fprintf(f, '%s\n',strcat('<track id="',num2str(counterTrack),'">'));
    for counterHop=1:handles.distanceNetwork.numHops(counterTrack)
        % loop to find the positions of each element of the track
        currentObject   = handles.nodeNetwork(handles.finalNetwork(counterHop,counterTrack),:);
        fprintf(f, '%s\n',strcat(   str1,num2str(currentObject(5)-1),...   % time
                                    str2,num2str(currentObject(1)),...   % X
                                    str3,num2str(currentObject(2)),...   % Y
                                    str4,num2str(currentObject(3)),...   % Z
                                    str5));
    end
    % end the track add labels
    fprintf(f, '%s\n','</track>');
end

% print the closing IDs and close file
fprintf(f, '%s\n','</trackgroup>');
fprintf(f, '%s\n','<linklist/>');
fprintf(f, '%s\n','</root>');  

fclose(f);



