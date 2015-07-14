function [tailBoundary,dataEdge,dataEdgesSum] =findTailBoundary (dataIn)
%function [tailBoundary,dataEdge,dataEdgesSum] =findTailBoundary (dataIn)
%
%-----------------------------------------------------------------
% findTailBoundary  fund tail bundary
%
%       INPUT
%         dataIn:           The greyscale data of a confocal scan of a
%                           zebra fish, separate of the fluorescence may be
%                           r x c x l x time.
%       OUTPUT
%         tailBoundary:     the tal boundary with same dimensions
%         dataEdge:
%         dataEdgesSum:
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


if ~exist('dataIn','var'); tailBoundary=[]; disp('no input data'); return; end
%%
%-----dimensions will determine if it is a 2D, 3D, 4D or a string that specifies a folder
if isa(dataIn,'char')
    folderName                                  = dataIn;
    if ~strcmp(folderName(end),'/');
        folderName                              = strcat(folderName,'/');
    end
    if isdir(folderName)
        dir1                                    = dir(strcat(folderName,'*.mat'));
        numImages                               = length(dir1);
        
        if (numImages>0)
            %%
            stepImages                              = floor(numImages/50); %#ok<NASGU>
            folderName                                  = dataIn;

            dir1                                    = dir(strcat(folderName,'*.mat'));
            numImages                               = length(dir1);
            stepImages                              = floor(numImages/50);

            %scan every 100 frames to create a edge image
            tempdata                            = load(strcat(folderName,dir1(1).name));
            [rows,cols,levs,timeFrames]         = size(tempdata.dataR); %#ok<NASGU>
            dataEdge(rows,cols,numel(1:stepImages:numImages)) = 0;
            k=1;
            %This counter will scan the data, skipping every 50 to obtain the edges at each step.
            %the idea is that strong edges will repeat and others will not and thus later a radon will
            %find the consistent edges to create a rough boundary and discard everyghin outside the tail
            for counterImages=1:stepImages:numImages
                tempdata                        = load(strcat(folderName,dir1(counterImages).name));
                dataEdge(:,:,k)                 = edge(tempdata.dataR(:,:,1+floor(levs/2)),'canny',[],1.5);
                k=k+1;
            end
            

            %%
            % Radon transform will find strongest lines, these are bound to be in the centre of the
            % notochord, then that will get the central orientation, get two areas of interest, one
            % oriented with this main orientation and another one at 90 degrees for the edge of the tail
            dataEdgesSum                            = sum(dataEdge(:,:,:),3);
            [dataRadon0]                            = radon(dataEdgesSum);
            [rowsRad,colsRad]                           = size(dataRadon0);
            %
            dataToClear                                 =zeros(rowsRad,colsRad);
            dataToClear(floor(0.33*rowsRad):ceil(0.66*rowsRad),:)      = 1;
            %
            peaksMain                               = houghpeaks(round(dataRadon0.*dataToClear),8,'NHoodSize',5*[1  1]);
            peaksMainL                              = min(peaksMain(:,2));
            peaksMainH                              = max(peaksMain(:,2));
            peaksMainA                              = mean(peaksMain(:,2));

            % get main orientation range and same for tail
            if peaksMainA>90; rotAngle =-90; else rotAngle = +90; end
            regionAngleToKeep1=(  peaksMainL-15:peaksMainH+15); %#ok<NASGU>
            regionAngleToKeep2= rotAngle+(peaksMainL-15:peaksMainH+15); %#ok<NASGU>

            %rotate in case the main angle is close to 0 or 180
            regionAngleToKeep2(regionAngleToKeep2<1) = regionAngleToKeep2(regionAngleToKeep2<1)+180;
            regionAngleToKeep2(regionAngleToKeep2>180) = regionAngleToKeep2(regionAngleToKeep2>180)-180;

            %crop any angles that escape [1-180]
            regionAngleToKeep1(regionAngleToKeep1>180)  = [];
            regionAngleToKeep2(regionAngleToKeep2>180)  = [];
            regionAngleToKeep1(regionAngleToKeep1<1)    = [];
            regionAngleToKeep2(regionAngleToKeep2<1)    = [];

            %

            dataToClear                                 =zeros(rowsRad-1,colsRad);
            dataToClear(1:floor(0.95*min(peaksMain(:,1))),regionAngleToKeep1)      = 1;
            dataDiffRadon                               = round(abs(diff(imfilter(dataRadon0,gaussF(5,5,1),'replicate'))));
            dataInAngles1                               =dataDiffRadon.*dataToClear;

            dataToClear                                 =zeros(rowsRad-1,colsRad);
            dataToClear(ceil(1.05*max(peaksMain(:,1))):end,regionAngleToKeep1)      = 1;
            
            dataInAngles2                               =dataDiffRadon.*dataToClear;
            
            dataInAngles3                            = dataDiffRadon;
            dataInAngles3(:,setxor((1:180),regionAngleToKeep2))      = 0;

            % now find the peaks, there should be 3 two for the sides that run parallel to the notochord and one
            % for the edge, perpendicular to the notochord, create a dummy variable to have for iradon
            %%
            peaks1                                   = houghpeaks(dataInAngles1,5,'NHoodSize',5*[1  1]);
            peaks2                                   = houghpeaks(dataInAngles2,5,'NHoodSize',5*[1  1]);
            peaks3                                   = houghpeaks(dataInAngles3,5,'NHoodSize',5*[1  1]);
            [minPeak1,indPeak1]                      = min(peaks1(:,1));
            [maxPeak2,indPeak2]                      = max(peaks2(:,1));
            [maxPeak3,indPeak3]                      = max(abs(peaks3(:,1)-(rowsRad/2)));

            peaks0=[peaks1(indPeak1,:);peaks2(indPeak2,:);peaks3(indPeak3,:)];
            %%
            dataPeaks                               =zeros(size(dataInAngles1));

            for counterPeaks = [1 2 3]
                dataPeaks(peaks0(counterPeaks,1),peaks0(counterPeaks,2))=1;
            end
            %obtain the inverse transform of the peaks0 and dilate to split into regions
            dataPeaksI                              = imdilate(iradon(dataPeaks,(0:179))>0,[1 1;1 1]);
            dataPeaksI2                             = bwlabel(dataPeaksI(2:end-1,2:end-1)==0);
            dataPeaksPr                             = regionprops(dataPeaksI2,'area');

            %%
            [maxA,indMaxA]                          = max([dataPeaksPr.Area]);
            dataPeakCentral                         = (dataPeaksI2==indMaxA);
            SE_Circ                                 = strel('disk',9);
            dataPeakCentralDil                      = imdilate(dataPeakCentral,SE_Circ);
            
            %%
            figure
            surfdat(dataPeakCentralDil.*dataEdgesSum)
        end
    end
end
%%
tailBoundary = dataPeakCentralDil;
