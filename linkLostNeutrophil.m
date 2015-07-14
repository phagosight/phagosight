function [numObjectsChanged,dataL,handles] = linkLostNeutrophil(handles,dataRe)
%function [numObjectsChanged,dataL,handles] = linkLostNeutrophil(handles,dataRe)
%
%--------------------------------------------------------------------------
% linkLostNeutrophil links cells that may be lost due to changes of intensity  
%      will have to loop all the dataRs, SAVE on the previously created dataL
%      folders and change the handles
%
%       INPUT
%         handles:              handles struct including nodeNetwork
%         dataRe:               path to reduced data
%
%       OUTPUT
%         numObjectsChanged:	number of linked neutrophils
%         dataL:                path to labelled data
%         handles:              updated handle struct
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
% 
%%
%% Find all the objects that disappear at one frame, then reappear after one black
recordLinks=[];
numChange = inf;
while numChange>0
    [handles,numChange,recordLinksT]=linkEndStartTracks(handles,4.3,2);
    recordLinks = [recordLinks;recordLinksT]; %#ok<AGROW>
end

%  Loops  over the FRAMES
numFramesToProcess                              = size(recordLinks,1);

% Generate the FILE NAME for the labelled folders
dataLa                                           = strcat(dataRe(1:end-2),'La');

% --------------------------- dataReName is not mat file, should be a folder with a)matlab files
dir2                                            = dir(strcat(dataLa,'/*.mat'));
%% Loop over the FRAMES
for counterDir=1:numFramesToProcess
    %%
    for counterT =[1 3 2]
        tempDirLa                                   = dir2(recordLinks(counterDir,1)+counterT-2+1,1).name;
        dataLaName                                  = strcat(dataLa,'/',tempDirLa);

        % now read the labelled data dataL for the times t-1, t, t+1
        dataFromFile                                =  load(dataLaName);
        if isfield(dataFromFile,'dataL')
            dataLa2                                 = dataFromFile.dataL;
        else
            namesF                                  = fieldnames(dataFromFile);
            dataLa2                                 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
        end
        dataLaT(:,:,:,counterT)                       = dataLa2; %#ok<AGROW>
    end
    
    [rows,cols,levs]                         = size(dataLa2); %#ok<NASGU>
    
    %detect limits of labelled data in the t-1 frame
    initC1                                   = max(1,-1+floor(handles.nodeNetwork(recordLinks(counterDir,4),15)));
    finC1                                    = min(cols,1+ceil(handles.nodeNetwork(recordLinks(counterDir,4),15)+handles.nodeNetwork(recordLinks(counterDir,4),18)));
    initR1                                   = max(1,-1+floor(handles.nodeNetwork(recordLinks(counterDir,4),16)));
    finR1                                    = min(rows,1+ceil(handles.nodeNetwork(recordLinks(counterDir,4),16)+handles.nodeNetwork(recordLinks(counterDir,4),19)));
    %detect limits of labelled data in the t+1 frame
    initC2                                   = max(1,-1+floor(handles.nodeNetwork(recordLinks(counterDir,5),15)));
    finC2                                    = min(cols,1+ceil(handles.nodeNetwork(recordLinks(counterDir,5),15)+handles.nodeNetwork(recordLinks(counterDir,4),18)));
    initR2                                   = max(1,-1+floor(handles.nodeNetwork(recordLinks(counterDir,5),16)));
    finR2                                    = min(rows,1+ceil(handles.nodeNetwork(recordLinks(counterDir,5),16)+handles.nodeNetwork(recordLinks(counterDir,4),19)));


    initR                                   = min(initR1,initR2);
    initC                                   = min(initC1,initC2);
    finR                                    = max(finR1,finR2);
    finC                                    = max(finC1,finC2);

    %this only works if there is NOTHING in the same region at time t
    if sum(sum(sum(dataLaT(initR1:finR1,initC1:finC1,:,2))))==0
        highestCurrentLabel                     = max(max(max(dataLaT(:,:,:,2))));

        dataLaT(initR:finR,initC:finC,:,2)      = (highestCurrentLabel +1) * ((dataLaT(initR:finR,initC:finC,:,1))&(dataLaT(initR:finR,initC:finC,:,3)));
        dataL = dataLaT(:,:,:,2); %#ok<NASGU>
        save(dataLaName,'dataL');
        
    end
end
if isempty(recordLinks)
    numObjectsChanged=0;
else
    numObjectsChanged = size(recordLinks,1);
end



