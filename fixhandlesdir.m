function [newhandles] = fixhandlesdir(handles)
% fixhandlesdir. Changes values for handles.dataRe and handles.dataLa to 
% fit the platform and (if on Windows) the location of the data.
% 
% USAGE:
%           [newhandles] = fixhandlesdir(handles)
%

newhandles = handles;
splitla = strsplit(handles.dataLa, filesep);
splitre = strsplit(handles.dataRe, filesep);

if isdir(handles.dataLa)
    fprintf('%s: Folder names appear to be consistent. No changes done.\n',...
        mfilename);
    return
end

switch chooseplatform
    case 'win'
        [~,a] = system('wmic logicaldisk get caption');
        b = strsplit(a,':');
        b = b{1};
        A = strsplit(a,b(end-1));
        
        ix = 1;
        for i=1:length(A)
            if ~isempty(strfind(A{i},':'))
                dirName{ix} = strcat(A{i},'\');
                ix=ix+1;
            end
        end
        
        % test now if change is from windows to windows, or
        % something else to windows.
        if ~isempty(strfind(splitla{1}, ':'))
            % windows2windows
            for ix=1:length(dirName)
                testdirname = fullfile(dirName{ix}, ...
                    joindirname(splitla(2:end)));
                if isdir(testdirname)
                    newhandles.dataLa = testdirname;
                    newhandles.dataRe = fullfile(dirName{ix},....
                        joindirname(splitre(2:end)));
                    break;
                end
            end
        else
            newdirLa = uigetdir('.', 'Select dataLa folder');
            newdirRe = uigetdir(newdirLa, 'Select dataRe folder');
            
            newhandles.dataLa = newdirLa;
            newhandles.dataRe = newdirRe;
        end
    case 'linux'
        newdirLa = uigetdir('.', 'Select dataLa folder');
        newdirRe = uigetdir(newdirLa, 'Select dataRe folder');
        
        newhandles.dataLa = newdirLa;
        newhandles.dataRe = newdirRe;
        
    case 'mac'
        newdirLa = uigetdir('.', 'Select dataLa folder');
        newdirRe = uigetdir(newdirLa, 'Select dataRe folder');
        
        newhandles.dataLa = newdirLa;
        newhandles.dataRe = newdirRe;
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