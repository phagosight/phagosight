function [dataR,handles] = reduceNeutrophil(dataIn,numLevels,handles)
%function [dataR,handles] = reduceNeutrophil(dataIn,numLevels,handles)
%
%--------------------------------------------------------------------------
% reduceNeutrophil  reduces spatial resolution of data
%   
%       INPUT
%         dataIn:           image to be reduced and labelled or path to the
%                           folder containing data to reduce.
%
%         numLevels:        levels to reduce in a pyramid, i.e. level = 1, 
%                           will reduce 1000 x 1000 x 5 to 500 x 500 x 5, 
%                           level =2 to 250 x 250 x 5, etc. This reduces 
%                           noise and computational complexity at the 
%                           expense of spatial resolution.
%
%         handles:          handles struct
%
%       OUTPUT
%         dataR:            data reduced in dimensions by uniform pyramid
%
%         handles:          updated handles struct           
%          
%--------------------------------------------------------------------------
%
%     Copyright (C) 2012  Constantino Carlos Reyes-Aldasoro
%
%     This file is part of the PhagoSight package.
%
%     The PhagoSight package is free software: you can redistribute it and/
%     or modify it under the terms of the GNU General Public License as 
%     published by the Free Software Foundation, version 3 of the License.
%
%     The PhagoSight package is distributed in the hope that it will be 
%     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with the PhagoSight package.  If not, see:
%
%                   <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
%
% This m-file is part of the PhagoSight package used to analyse fluorescent 
% phagocytes as observed through confocal or multiphoton microscopes.  For 
% a comprehensive user manual, please visit:
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
% accuracy, completeness, or usefulness of any information, or method in  
% the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------
if ~exist('numLevels','var')
    numLevels=2;
end

% Reduce complexity of data set,
if isa(dataIn,'char')
    % dataInName is not mat file, should be a folder with a)matlab files
    %disp('Read matlab files from folder, reduce and save in a new folder')

    dir1 = dir(strcat(dataIn,'/*.mat'));

    dataOutName = strcat(dataIn(1:end-2),'Re');
    mkdir(dataOutName)

    numFrames = size(dir1,1);
    for counterDir=1:numFrames
      tempDir = dir1(counterDir).name;
      dataInName = strcat(dataIn,'/',tempDir);
      dataOutName1 = strcat(dataOutName,'/',tempDir);
      
      dataFromFile = load(dataInName);
      if isfield(dataFromFile,'dataIn')
        dataIn2 = dataFromFile.dataIn;
      else
        namesF = fieldnames(dataFromFile);
        dataIn2 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
      end
      dataR = reduceNeutrophil(dataIn2,numLevels);
      save(dataOutName1,'dataR');
    end
    
    dataR = dataOutName;
else
    if numLevels==0
      dataR = (imfilter(double(dataIn),ones(3)/9,'replicate'));
      %dataR = (imfilter(double(dataIn),fspecial('gaussian'),'replicate'));
    else
      dataR = reduceu(dataIn,numLevels);
    end
end

% Regular dimension check and definition of time frames
if ~isa(dataR,'char')
  [handles.rows,handles.cols,handles.levs] = size(dataR);
end
