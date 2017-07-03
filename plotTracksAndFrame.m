function plotTracksAndFrame(handles, whichFrame, winsize, options)
% plotTracksAndFrame. Plots the with plotTracks function around a
% frame(whichFrame) and zoomed in to the time axis (z-axis) in a window
% determined by (winsize). The options can be set to determine the
% parameter typeOfPlot of function plotTracks, as well as the minimum
% number of hops (numHops) necessary to plot the track or the type of data
% that will be displayed ('La' or 'Re').
%
% USAGE:
%       plotTracksAndFrame(handles, whichFrame)
%       plotTracksAndFrame(handles, whichFrame, winsize)
%       plotTracksAndFrame(handles, whichFrame, winsize, options)
%
% INPUT:
%               handles := structure containing
%                   nodeNetwork : [numRBC detected x 12 params]
%                   finalNetwork: either 1 track or [depth of tracks x numTracks]
%                   dataRe      : string with path to Reduced data (mat_Re)
%                   dataLa      : string with path to Labelled data (mat_La)
%               whichFrame := Frame between 1 and handles.numFrames to
%                           display on top of the tracks.
%
%               winsize := size of the window that the tracks will zoom
%                           into the time (z) axis.
%
%               options :=  structure including at leaast one of the next
%                           fields:
%                   numHops    : minimum number of "hops" from
%                       handles.distanceNetwork.numHops to plot any track.
%                   typeOfPlot : number from 1 - 14 for PLOTTRACKS function
%                       or 'none' to plot tracks in green when
%                       handles.distanceNetwork.numHops > numHops
%                       or red when
%                       handles.distanceNetwork.numHops <= numHops
%                   typeOfData : 'La', 'Re' or 'none' to determine which
%                       frame will be placed alongside the tracks.
%
% see also plotTracks
%

if ~isdir(handles.dataLa)
    %fprintf('%s: Checking conistency in dataRe/dataLa folder paths.\n', mfilename);
    handles = fixhandlesdir(handles);
end

if nargin<4
    typeOfPlot = [];
    typeOfData = 'none';
    numHops = 50;
    tracks2plot = 'all';
elseif nargin<3
    winsize = 10;
    [typeOfPlot, typeOfData, numHops, tracks2plot] = getOptions(options);
else
    [typeOfPlot, typeOfData, numHops, tracks2plot] = getOptions(options);
end

switch typeOfData
    case 'La'
        bnames = dir(handles.dataLa);
        bnames(1:2) = [];
        bnames = {bnames.name};
        currData = load(fullfile(handles.dataLa,bnames{whichFrame}));

        X = getfield(currData, 'dataL');

    case 'Re'
        bnames = dir(fullfile(handles.dataRe, '*.mat'));
        bnames = {bnames.name};
        currData = load(fullfile(handles.dataRe,bnames{whichFrame}));
        structnames = fieldnames(currData);

        X = getfield(currData, structnames{1});

    case 'none'

    otherwise
        fprintf('%s: ERROR non recognised options.typeOfData: \n\t%s\n',...
            mfilename, typeOfData);
        fprintf('%s: try: "Re", "La".\n', mfilename);
        X = zeros(handles.rows, handles.cols);

end

[xx,yy] = meshgrid(1:handles.cols, 1:handles.rows);
zz = whichFrame.*ones(size(xx));

% z-axis values
zmin = max(whichFrame - winsize, 1);
zmax = min(whichFrame + winsize, handles.numFrames);

cla;
if isempty(typeOfPlot)
    titlestr = strcat('Short and long tracks. Centre frame =',...
        32, num2str(whichFrame), 32, 'Windows size =', 32, ...
        num2str(winsize));
    title(titlestr);
    plotTracks(handles, 11, find(handles.distanceNetwork.numHops>numHops));
    hold on;
    plotTracks(handles, 12, find(handles.distanceNetwork.numHops<=numHops));
    if ~strcmp(typeOfData, 'none')
        surface(xx,yy,zz,X, 'FaceColor', 'texturemap', ...
            'EdgeColor','none','CDataMapping','direct');
        alpha(0.75);
    end
    hold off;
    axis([1 handles.cols 1 handles.rows zmin zmax]);
else
    titlestr = strcat('Tracks plot. Centre frame =',...
        32, num2str(whichFrame), 32, 'Windows size =', 32, ...
        num2str(winsize));
    title(titlestr);
    plotTracks(handles, typeOfPlot, ...
        find(handles.distanceNetwork.numHops>numHops));
    if ~strcmp(typeOfData, 'none')
        hold on;
        surface(xx,yy,zz,X, 'FaceColor', 'texturemap', ...
            'EdgeColor','none','CDataMapping','direct');
        alpha(0.75);
        hold off;
    end
    axis([1 handles.cols 1 handles.rows zmin zmax]);
end

end

function [typeOfPlot, typeOfData, numHops, tracks2plot] = getOptions(s)
% Get options for this function
%
typeOfPlot = [];
typeOfData = 'La';
numHops = 50;
tracks2plot = 'all';

fnames = fieldnames(s);
for ix=1:length(fnames)
    name = fnames{ix};
    switch name
        case 'typeOfPlot'
            typeOfPlot = getfield(s, name);
            if ischar(typeOfPlot)
                typeOfPlot = str2num(typeOfPlot);
            end
        case 'typeOfData'
            typeOfData = getfield(s, name);
        case 'numHops'
            numHops = getfield(s, name);
            if ischar(numHops)
                if strcmp(numHops, 'none')
                    numHops = 0;
                else
                    numHops = str2num(numHops);
                end
            end
        case 'tracks2plot'
            tracks2plot = s.(name);
            numHops = 1; 

        otherwise
            fprintf('%s: ERROR, incorrect option selected: %s is NOT defined\n',...
                mfilename, upper(name));
    end
end
end
