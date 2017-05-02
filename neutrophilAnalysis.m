function [handles]= neutrophilAnalysis(dataInName,numLevelsToReduce,...
                                        dataStructure,thresLevels,minBlob)
%[handles]= neutrophilAnalysis()
%[handles]= neutrophilAnalysis(dataInName)
%[handles]= neutrophilAnalysis(dataInName,numLevelsToReduce)
%[handles]= neutrophilAnalysis(dataInName,numLevelsToReduce,...
%                               dataStructure)
%[handles]= neutrophilAnalysis(dataInName,numLevelsToReduce,...
%                               dataStructure,thresLevels)
%[handles]= neutrophilAnalysis(dataInName,numLevelsToReduce,...
%                               dataStructure,thresLevels,minBlob)
%
%
%--------------------------------------------------------------------------
% neutrophilAnalysis is the main routine for the  neutrophil analysis
%     All other subroutines are called from here. It computes the handles
%     struct containing information of the tracks detected in the input.
%
%       INPUT
%         dataInName:           (a) the name of a folder with matlab files, 
%                                   one per time frame, each matrix is 3D
%                               (b) the name of a folder with folders with  
%                                   tiff files, one per frame.
%
%         numLevelsToReduce:    for tiff files or not reduced mat files, id
%                               indicates the number of levels that the data 
%                               will be reduced in a pyramid, i.e. level = 1, 
%                               will reduce 1000 x 1000 x 5 to 500 x 500 x 5, 
%                               level =2 to 250 x 250 x 5, etc. This reduces 
%                               noise and computational complexity at the 
%                               expense of spatial resolution.
%                               Level = 0 leaves the data unchanged. 
%                               Default is 0.
%
%         dataStructure:        6x1 vector containing the slice's first and
%                               last indexes of the Green, DIC and Red Channel
%                               Example: [6 10 1 5 11 15]' means slices 1 to 5
%                               correspond to the DIC channel, slices 6 to 10
%                               to the Green channel and, finally, slices from
%                               11 to 15 to the Red channel.
%
%         thresLevels:          1x2 vector containing minimum and maximum
%                               values to threshold data.
%
%         minBlob:              Blobs smallest than this will be considered 
%                               noise. 
%
%
%       OUTPUT
%         handles:              a struct with fields that described the data sets 
%                               in terms of the tracking of neutrophils, and 
%                               many characteristics. The handles are further used 
%                               to extract statistics, plot the tracks, etc.
%
%                                  |  1. X |
%                                  |  2. Y |- Position of centroid
%                                  |  3. Z |
%                                  |  4. Distance to closest Neutrophil.
%                                  |  5. Time Slice
%                                  |  6. ID
%                                  |  7. Parent (0 if initial node of track)
%                                  |  8. Child (0 if no children)
%                                  |  9. Velocity
%             handles.nodeNetwork -| 10. Volume
%                                  | 11. Label at t
%                                  | 12. Keyhole Region (1-circle, 2-wedge)
%                                  | 13. Track (identifier)
%                                  | 14. Final Label
%                                  | 15. X init |
%                                  | 16. Y init |
%                                  | 17. Z init |
%                                  | 18. size X |- Bounding box
%                                  | 19. size Y |
%                                  | 20. size Z |
%
%             handles.finalNetwork [ one column for every track, ...
%                                      each column contains the nodes ...
%                                      that belong to that track]
%
%                               The processing also generates many folders:
%                               (a) dataInName_mat_Or The original data (no 
%                                   reduction) but in matlab format
%                               (b) dataInName_mat_Re Data reduced
%                               (c) dataInName_mat_La Data as labelled images 
%                                   after segmentation
%                               (d) dataInName_mat_Ha Data handles and other 
%                                   final data (wound regions, etc.) 
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
%     along with the PhagoSight package.  If not, see 
%
%            <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
%
% This m-file is part of the PhagoSight package used to analyse
% fluorescent phagocytes as observed through confocal or
% multiphoton microscopes.  For a comprehensive manual, please visit:
%
%            <http://www.phagosight.org.uk>
%
% Please feel welcome to use, adapt or modify the files. If you can improve
% the performance of any other algorithm please contact us so that we can
% update the package accordingly.
%
%--------------------------------------------------------------------------
%
% The authors shall not be liable for any errors or responsibility for the 
% accuracy, completeness, or usefulness of any information, or method in 
% the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

%test the path to the PhagoSight package it must be able to detect but the 
% m-files of PhagoSight. 
% If they are not readable the path must be set first. 
testPhagoSight  = which('neutrophilAnalysis');

if (isempty(testPhagoSight))
    %the Path has not been set properly, display legend and exit.
    disp('')
    disp('---------------------------------------------------')
    disp('   PhagoSight has not been properly configured')
    disp('---------------------------------------------------')
    disp('PhagoSight requires a "Path" to all the relevant ')
    disp('files. To do this, go to: ')
    disp('    FILE > PATH > ADD FOLDER WITH SUBFOLDERS')
    disp('And then select the folder "PhagoSight". ')
    disp('Matlab will add a number of folders to the list ')
    disp('in the window. Click SAVE and CLOSE.')
    disp('After this PhagoSight is ready to be used.')
    disp('You only need to do this once.')
    disp('---------------------------------------------------')
    disp('')   
    handles = [];
    return;
else
    clear testPhagoSight
end 

% Parsing of input arguments, if no data is received, call readNeutrophils
% that will open a new window to locate the folder where the data is stored
if (nargin == 0)
    %No parameters, readNeutrophil shows a window when user can select
    %the data folder.
    [dataIn,handles] = readNeutrophils();
    if isempty(handles)
            return;
    else
        STR = [num2str(handles.numFrames) ' frames detected.'];
        disp(STR);
    end
elseif (isempty(dataInName))
    disp('dataInName is empty, calling readNeutrophil to select data source...');
    [dataIn handles] = readNeutrophils();
    if isemtpy(handles)
        return;
    else        
        STR = [num2str(handles.numFrames) ...
               ' frames detected.'];
        disp(STR);
    end
else
    %Al least 1 parameter. The first parameter must be a folder
    %containing either mat files or folders with tiff files; for the
    %former, dataInName is the path to a folder named *_mat_*.
    if (strfind(dataInName,'_mat_'))
        %Folder containing mat files, determine how many files/frames.
        %As information about the number of levels to reduce is not
        %provided, a default value 1 is used.
        lstFiles = dir(strcat(dataInName,'/*.mat'));
        if isempty(lstFiles )
            disp('Folder not found');
            handles = [];
            return;
        else
            %                 nbFolders                           = 0;
            %                 while(lstFiles(nbFolders+1).isdir)
            %                     nbFolders = nbFolders + 1;
            %                 end
            %                 handles.numFrames = size(lstFiles,1)-nbFolders;
            handles.numFrames = size(lstFiles,1);
            dataIn = dataInName;
        end
    else
        %Folder contains a tree of folders, one per frame,
        %containing tiff files, each one representing a plane of the
        %data volume. readNeutrophil allow to read those files and
        %parse them to matlab data files (mat file)
        [dataIn,handles] = readNeutrophils(dataInName); 
        
        if isempty(dataIn)
            disp('Folder not found');
            handles = [];
            return;
        else
            STR = [num2str(handles.numFrames) ...
                   ' frames detected.'];
            disp(STR);
        end
    end
end

%%Verifying data was acquired.
if (isempty(dataIn))
    disp('Error determining source of data.');
    handles = [];
    return;
end

%Reduce computational complexity by averaging neighbouring pixels
%in a uniform pyramid, 2 levels would reduce from 1000x1000 to 250x250
%reduceNeutrophils will read the folder and then iteratively reduce
%size in each frame and save if dataIn is a pathname, or will
%return the data as a matrix
if (~exist('numLevelsToReduce','var'))
    numLevelsToReduce = 0;
else
    if ischar(numLevelsToReduce)
        numLevelsToReduce = str2num(numLevelsToReduce);
    end
end
if exist('minBlob','var')
    handles.minBlob = minBlob;
end

if (isa(dataIn,'char'))
    %Removing separator character from path
    if (strcmp(dataIn(end),'/')) 
        dataIn = dataIn(1:end-1); 
    end
    %The processingLevel = 1 was eliminated 'cause that processing task
    %was already executed in the parsing some lines above
    switch (dataIn(end-6:end))
      case '_mat_Or'
        processingLevel                 = 2;
      case '_mat_Re'
        processingLevel                 = 3;
      case '_mat_La'
        processingLevel                 = 4;
      case '_mat_Ha'
        processingLevel                 = 5;
      otherwise
        processingLevel                 = 1;
    end
    
    dataRe = strcat(dataIn(1:end-2),'Re');
    dataLa = strcat(dataIn(1:end-2),'La');
        
    if (processingLevel < 3)
        %ProcessingLevel = 2
        % fileName_mat_Or
        % dataIn contains original matlab data: read and reduce
        disp(['Read matlab files from folder,',32,...
              'reduce and save in a new folder *_mat_Re']);
        %if (numLevelsToReduce == 0)
        %    dataR = dataIn; 
        %else
        
        [dataR,handles] = reduceNeutrophil(dataIn,...
                                           numLevelsToReduce,...
                                           handles);
        %end
    end
    
    if (processingLevel < 5)
        %ProcessingLevel = 3
        %ProcessingLevel = 4
        
        if (~exist('handles','var') || ...
            (exist('handles','var') && ~isfield(handles,'ChannelDistribution')))
            % Determine the DATA STRUCTURE, it may be that the data is:
            %       Top 50% fluorescent - bottom 50% DIC
            %       Top 50% fluorescent - bottom 50% fluorescent
            %       Top X%  fluorescent - bottom Y%  DIC
            %       Only fluorescent 
            %  The distribution of the channels will be contained in the new field:
            %      handles.ChannelDistribution  [Green_init
            %      Green_fin   DIC_init DIC_fin Red_init Red_fin]
            if exist('dataStructure','var')
                %In case the channel distribution was provided by the
                %user as an argument for the function
                
                [handles,dataRe2] = determineDataStructure(dataRe,dataStructure);
            else
                %Ask the user to provide the channel distribution
                
                [handles,dataRe2] = determineDataStructure(dataRe);
            end
            if (~isfield(handles,'ChannelDistribution'))
                disp('An ERROR occured while determining the Channel Distribution');
                handles=[];
                return;
            end
        end
        %%Data structure not present at this stage indicates that
        %%there was a cancel or an error in the process, exit
        if (isempty(handles.ChannelDistribution))
            disp('Channel Distribution not assigned.');
            handles = [];
            return;
        end
        if exist('minBlob','var')
            handles.minBlob = minBlob;
        end
        
        
        %Saving location of reduced and labelled data into the
        %handles variable.
        handles.dataRe = dataRe;
        handles.dataLa = dataLa;
            
        if (~isfield(handles,'thresLevel'))
            % ---------- threshold and label the neutrophils. ------------------
            % To increase the accuracy for the segmentation a
            % double thresholding  is performed, small blobs are removed as noise.
            % Thresholds are automatically calculated with Otsu
            % and then manually verified in setNeutrophilHandes
            
            if exist('thresLevels','var')
                %User provided the threshold levels as an argument of
                %the function
                handles = setNeutrophilHandles(dataRe2,handles,thresLevels);
            elseif processingLevel==4
                handles = setNeutrophilHandles(dataRe2,handles,[], false);
            else
                %Calculate and validate threshold levels
                % TO REMOVE USER INPUT! USERINPUT
                handles = setNeutrophilHandles(dataRe2,handles,[]);
            end
        end
    end
    
    if (processingLevel < 4)
        % fileName_mat_Re
        [dataL,handles] = thresNeutrophil(dataRe,handles);
        
        % -------- Split objects that are outliers by size  ----------------------
        [handles,numObjectsToSplit] = splitLargeNeutrophil(handles,dataLa,[]);
    end
    if (processingLevel < 5)
        
        if (~isfield(handles,'rows'))
            handles = setNeutrophilHandles(dataIn,handles,[],false);                
        end
        
        if (~isfield(handles,'nodeNetwork') || size(handles.nodeNetwork,1) > 0)
            
            % fileName_mat_La
            % dataIn contains labelled data: read and track
            if ~exist('dataL','var');
                dataL = dataIn;
            end
            disp('Initial Tracking process')    
            handles = trackNeutrophil(dataL,handles,0);
            
            %-- Find cases where a neutrop disappears 1 frame, redraw and re-link --
            numObjectsChanged = linkLostNeutrophil(handles,dataRe);
            if numObjectsChanged>0
                handles2 = rmfield(handles,{'nodeNetwork','neighNetwork',...
                                    'linkedNetwork','finalNetwork',...
                                    'distanceNetwork','finalLabel',...
                                    'reducedNetwork'});
                %-- with the new "lost" neutrops, track again -------------
                handles2 = trackNeutrophil(dataL,handles2,0);
            else
                handles2 = handles;
            end
            %--  Split objects that merge by collisions ----------------------
            numObjectsSplit = splitNeutrophilCollision(handles2,dataRe);
            if numObjectsSplit>0
                handles3 = rmfield(handles2,{'nodeNetwork','neighNetwork',...
                                    'linkedNetwork','finalNetwork',...
                                    'distanceNetwork',...
                                    'finalLabel','reducedNetwork'});
                %-- with the new "lost" neutrops, track for a final time  ------
                disp('Final Tracking process');
                handles3 = trackNeutrophil(dataL,handles3,0);
            else
                handles3 = handles;
            end
            %-- link tracks that are broken between consecutive frames --------

            handlesPreSplit = handles;
            handles = handles3;
            %save the handles as they will be used later on
            dataHa = strcat(dataIn(1:end-2),'Ha');
            mkdir(dataHa)
            save(strcat(dataHa,'/handles.mat'),'handles','handlesPreSplit');
            
        else
            %There are no nodes to process
            disp(...
                strcat('It was NOT possible to detect any node in the data,',...
                       'verify the threshold levels!'));
            handles = [];
            return;
        end
    end
    if processingLevel<6
        % fileName_mat_Ha
        % dataIn contains tracked data, only get the handles
        dataHa = strcat(dataIn(1:end-2),'Ha/handles.mat');
        load(dataHa);
        
    end
    dataHa = strcat(dataIn(1:end-2),'Ha/handles.mat');
    load(dataHa);
    dataR = strcat(dataIn(1:end-2),'Re');
    
    handles.numLevelsToReduce = numLevelsToReduce;
    % with the use of a wound region, the effective distances
    % are calculated, and with that the networks
    % of translation and exploration
    woundFile = strcat(dataIn(1:end-2),'Wound/woundRegion.mat');
    if (exist(woundFile,'file'))
        load(woundFile);
        handles = effectiveDistance(handles,woundR);
        dataF = strcat(dataIn(1:end-2),'Ha/finalHandles.mat');
        save (dataF,'handles','handlesPreSplit','woundRegion','dataR');       
    else
        disp('(Wound file NOT found.)');
    end
    disp('Processing has finished without errors.');
else
    disp('An ERROR has occur while processing the data.');
end
