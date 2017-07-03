function [newhandles] = fixhandlesdir(handles)
% fixhandlesdir. Changes values for handles.dataRe and handles.dataLa to
% fit the platform and (if on Windows) the location of the data.
%
% USAGE:
%           [newhandles] = fixhandlesdir(handles)
%

newhandles = handles;

fnames = fieldnames(handles);
fnames = fnames(contains(fnames,'data'));

if isdir(handles.(fnames{1}))
    fprintf('%s: Folder names appear to be consistent. No changes done.\n',...
        mfilename);
    return
end

switch chooseplatform
    case 'win'
        [~,volumes] = system('wmic logicaldisk get caption');
        b = strsplit(volumes,':');
        b = b{1};
        A = strsplit(volumes,b(end-1));

        ix = 1;
        clear volumes;
        for i=1:length(A)
            if ~isempty(strfind(A{i},':'))
                volumes{ix} = strcat(A{i},'\');
                ix=ix+1;
            end
        end
    case 'linux'
        volumes = dir('/home');
        volumes(1:2) = [];
        volumes = {volumes.name};
        
        for jx=1:length(volumes)
            volumes{jx} = fullfile('/media', volumes{jx});
        end
    case 'mac'
        volumes = dir('/Volumes');
        volumes(1:2) = [];
        volumes = {volumes.name};
        
        for jx=1:length(volumes)
            volumes{jx} = fullfile('/Volumes', volumes{jx});

        end
end

for jx=1:length(volumes)
    for kx=1:length(fnames)
        [~, splitpath] = comesfrom(handles.(fnames{kx}));
        dirtest = joindirname({volumes{jx}, splitpath});
        if isdir(dirtest)
            newhandles.(fnames{ix}) = dirtest;
            break;
        end
    end
end

end

function [platformused] = chooseplatform()
% return the platform in which phagosight is being run.
test = isunix + ismac;

switch test
    case 0
        platformused = 'win';
    case 1
        platformused = 'linux';
    case 2
        platformused = 'mac';
end
end

function [joineddir] = joindirname(splitdirname)
% join folder name from splitstr's output.
joineddir ='';
for jx=1:length(splitdirname)
    joineddir = strcat(joineddir, filesep, splitdirname{jx});
end
end

function [op, splitpath] = comesfrom(somedatapath)
% returns which platform does the
if contains(somedatapath, ':')
    op = 'win';
    if nargout > 1
        splitpath = strsplit(somedatapath, '\');
        splitpath(1) = [];
    end
else
    if contains(somedatapath, '/Volumes')
        op = 'mac';
        if nargout > 1
            splitpath = strsplit(somedatapath, '/');
            splidx = find(contains(splitpath, 'Volumes'));
            splitpath(1:splidx) = [];
        end
    else
        op = 'linux';
        if nargout > 1
            splitpath = strsplit(somedatapath, '/');
            splidx = find(contains(splitpath, 'media'));
            splitpath(1:(splidx+1)) = [];
        end
    end
end
end