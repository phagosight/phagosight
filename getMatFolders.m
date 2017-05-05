function [handlesDir] = getMatFolders(baseFileName)
%               GET *_MAT_* FOLDERS
%
% GETMATFOLDERS: Given a path to data (baseFileName), return a structure
% (handlesDir) that contains the paths to directories created by the
% phagosight package,
%
%  EXAMPLE: if baseFileName='/path/to/data', then the folders
% 'path/to/data_mat_Ha', 'path/to/data_mat_Or', 'path/to/data_mat_Re' and
% 'path/to/data_mat_La' are searched for and returnet in structure
% handlesDir.
%
% USAGE:
%           [handlesDir] = getMatFolders();
%           [handlesDir] = getMatFolders(baseFileName);
% INPUT
%               baseFileName := (string) Path to directory that contains
%               the data.
% OUTPUT
%               handlesDir   := (struct) Structure with the following
%               parts:
%                 handlesDir.pathtodir := pathtodir
%                 handlesDir.data := baseFileName
%                 handlesDir.dataOr := contains _mat_Or folder
%                 handlesDir.dataRe := contains _mat_Re folder
%                 handlesDir.dataLa := contains _mat_La folder
%

if nargin < 1 || ~isdir(baseFileName)
    baseFileName = uigetdir('.', 'Select folder with RAW data.');
    if ~isdir(baseFileName)
        fprintf('%s: Incorrect baseFileName variable. Not a directory.',...
            mfilename);
    end
elseif strcmp(baseFileName(end), filesep)
    baseFileName(end) = [];
end

handlesDir.pathtodir = [];
handlesDir.dataHa = [];
handlesDir.dataOr = [];
handlesDir.dataRe = [];
handlesDir.dataLa = [];

if contains(baseFileName, '_mat_')
    indx = strfind(baseFileName, '_mat_');
    baseFileName = baseFileName(1:indx-1);
    lsdir = dir(strcat(baseFileName(1:indx-1), '_mat_*'));
else
    lsdir = dir(strcat(baseFileName, '_mat_*'));
end

if ~isempty(lsdir)
    indx2path = strfind(baseFileName, filesep);
    indx2path = indx2path(end);
    pathtodir = baseFileName(1:indx2path);
    handlesDir.pathtodir = pathtodir;
    handlesDir.data = baseFileName(indx2path+1:end);
    foldernames = {lsdir.name};
    for ix=1:length(foldernames)
        switch foldernames{ix}(end-1:end)
            case 'Ha'
                handlesDir.dataHa = fullfile(foldernames{ix});
            case 'Or'
                handlesDir.dataOr = fullfile(foldernames{ix});
            case 'Re'
                handlesDir.dataRe = fullfile(foldernames{ix});
            case 'La'
                handlesDir.dataLa = fullfile(foldernames{ix});
        end
    end
end


