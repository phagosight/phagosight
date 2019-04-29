function [handles, prevhandles] = removeMultipleTracks(handles,tracks2rm,woundRegion)
% REMOVEMULTIPLETRACKS removes multiple tracks from the handles structure. It
% uses the function removeOneTrack from the phagosight package.
% 
% USAGE:
%       [handles] = removeMultipleTracks(handles,tracks2rm)
%       [handles] = removeMultipleTracks(handles,tracks2rm, woundRegion)
%       [handles, prevhandles] = removeMultipleTracks(...)
% 
% INPUT:
%             handles := handles struct
%           tracks2rm := label (id numbers) of the tracks to remove
%         woundRegion := wound region matrix (to recalculate metrics)
%
% OUTPUT:
%             handles := updated handles
%         prevhandles := (optional) handles before processing
% 
% see also removeOneTrack
% 
if nargin<2 || any(tracks2rm<1)||length(tracks2rm)>size(handles.finalNetwork,2)
    fprintf('%s: Please verify that the input parameters are correct.\n\n',...
        mfilename);
    %help removeMultipleTracks;
    return;
end

if nargin < 3
    woundRegion = zeros(handles.rows,handles.cols);
end

if nargout > 1
    prevhandles = handles;
end

numtracks = length(tracks2rm);
for ix = numtracks:-1:1
    thistrack = tracks2rm(ix);
    try
        handles = removeOneTrack(handles, thistrack, woundRegion);
    catch e
        fprintf('%s: Error processing track %d. Moving on to next tracks.\n',...
            mfilename,tracks2rm(ix));
        disp(e.message);
    end
end


    


    
