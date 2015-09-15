function [dataL,handles,statsData]=thresNeutrophil(dataIn,handles)
%function [dataL,handles,statsData]=thresNeutrophil(dataIn,handles)
%
%--------------------------------------------------------------------------
% thresNeutrophil threshold data to generate labelled images containing only
%     neutrophils
%
%       INPUT
%         dataIn:       image to be thresholded and labelled
%         handles:      struct with all the parameters
%
%       OUTPUT
%         handles:      updated handles struct.
%         dataL:        data labelled from the threshold and small regions 
%                       removed
%         statsData:    statistics about labelled regions
%
%--------------------------------------------------------------------------
%
%     Copyright (C) 2012  Constantino Carlos Reyes-Aldasoro
%
%     This file is part of the PhagoSight package.
%
%     The PhagoSight package is free software: you can redistribute it 
%     and/or modify it under the terms of the GNU General Public License 
%     as published by the Free Software Foundation, version 3 of the License.
%
%     The PhagoSight package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with the PhagoSight package.  If not, see:
% 
%                   <http://www.gnu.org/licenses/>.
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
% accuracy, completeness, or usefulness of any information, or method in 
% the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------

% handles is not received, calculate the basic parameters
if ~exist('handles','var')
    maxData = max(dataIn(:));
    minData = min(dataIn(:));
    handles.thresLevel = round(0.5*maxData+0.5*minData);
    handles.thresLevel = 255*graythresh(uint8(dataIn(:)));
    handles.minBlob = 10;
    handles.ChannelDistribution = [1 size(dataIn,3) 0 0 0 0];
end

% two threshold levels are received in place of handles assign the threshold 
% levels and assume everything is a single channel in Green
if ~isstruct(handles)
    maxData = max(dataIn(:));
    minData = min(dataIn(:));
    tempValue = handles; clear handles;
    handles.thresLevel = [min(tempValue) max(tempValue)];
    handles.minBlob = 10;
    handles.ChannelDistribution = [1 size(dataIn,3) 0 0 0 0];
    
end


indexGreen = (handles.ChannelDistribution(1):handles.ChannelDistribution(2));
indexRed = (handles.ChannelDistribution(5):handles.ChannelDistribution(6));
indexDIC = (handles.ChannelDistribution(3):handles.ChannelDistribution(4));
    
if (indexRed~=0)&(numel(handles.minBlob)==1)
    handles.minBlob(2) = handles.minBlob(1);
end

if isa(dataIn,'char')
    % --------------------------- dataInName is not mat file, should be a folder with a)matlab files
    dir1 = dir(strcat(dataIn,'/*.mat'));
    dataOutName = strcat(dataIn(1:end-2),'La');
    mkdir(dataOutName)
    numFrames = size(dir1,1);
    for counterDir=1:numFrames
        tempDir = dir1(counterDir).name;
        dataInName =  strcat(dataIn,'/',tempDir);
        dataOutName1 =  strcat(dataOutName,'/',tempDir);
        
        dataFromFile =  load(dataInName);
        if isfield(dataFromFile,'dataIn')
            dataIn2 = dataFromFile.dataIn;
        else
            namesF = fieldnames(dataFromFile);
            dataIn2 = getfield(dataFromFile,namesF{1}); %#ok<GFLD>
        end
        %%        
        % handles.thresLevel can be a single value or a pair of values 
        % (high, low) if it is one, split into two to have separate 
        % thresholds
        if numel(handles.thresLevel)==1
            lowT_R = 0.95*(handles.thresLevel);
            highT_R = 1.05*(handles.thresLevel);            
            lowT_G = lowT_R;
            highT_G = highT_R;
        elseif numel(handles.thresLevel)==2
            lowT_R = min(handles.thresLevel);
            highT_R = max(handles.thresLevel);      
            lowT_G = lowT_R;
            highT_G = highT_R;
        elseif numel(handles.thresLevel)==4
            lowT_G = min(handles.thresLevel(1:2));
            highT_G = max(handles.thresLevel(1:2));            
            lowT_R = min(handles.thresLevel(3:4));
            highT_R = max(handles.thresLevel(3:4));                      
        else 
            handles=[]; return;
        end
        
        % handles.ChannelDistribution is used to determine which slices 
        % correspond to Green/Red/DIC
  
             if max(indexGreen)>0
                % Only process if there are GREEN channel slices
                 [GreenLowT,numGreenLowT] = ...
                     bwlabeln(dataIn2(:,:,indexGreen)>lowT_G);
                 
                 if numGreenLowT>0
                     GreenHighT = (dataIn2(:,:,indexGreen)>highT_G);
                     prodHighLowT = GreenLowT.*GreenHighT;
                     GreenLowT_2D = (sum(GreenLowT,3))>0;
                     clear rr cc
                     [rr,cc] = find(GreenLowT_2D);
                     rr = unique(rr);
                     cc = unique(cc);

                     neutropToKeep = unique(prodHighLowT(rr,cc,:));

                     [dataL_green,numNeutropGreen] = ...
                         bwlabeln(ismember(GreenLowT,neutropToKeep(2:end))); 
                 else
                     %This line assumes that there is only green channel
                     %dataL_green = zeros(size(dataIn2));
                     %Fixed here:
                     dataL_green = zeros([size(dataIn2,1) ...
                                            size(dataIn2,2) ...
                                            numel(indexGreen)]);
                     numNeutropGreen = 0;
                 end
             else
                 numNeutropGreen = 0;
             end

             if max(indexRed)>0
                 % Only process if there are RED channel slices
                 [RedLowT,numRedLowT] = ...
                     bwlabeln(dataIn2(:,:,indexRed)>lowT_R);
                 
                 if numRedLowT>0
                     RedHighT = (dataIn2(:,:,indexRed)>highT_R);
                     prodHighLowT = RedLowT.*RedHighT;
                     RedLowT_2D = (sum(RedLowT,3))>0;
                     clear rr cc
                     [rr,cc] = find(RedLowT_2D);
                     rr = unique(rr);
                     cc = unique(cc);
                     
                     neutropToKeep = unique(prodHighLowT(rr,cc,:));
                     
                     [dataL_Red,numNeutropRed] = ...
                         bwlabeln(ismember(RedLowT,neutropToKeep(2:end)));                   
                 else
                     %This line assumes that there is no green channel, 
                     %only red.
                     %dataL_Red = zeros(size(dataIn2));
                     %Fixed here:
                     dataL_Red = zeros([size(dataIn2,1) ...
                                        size(dataIn2,2) ...
                                        numel(indexRed)]);
                     numNeutropRed = 0;
                 end
             else
                 numNeutropRed = 0;
             end
            %   
            % now merge BACK into a single matrix if only one case
            if (max(indexDIC) > 0)
                dataL(:,:,indexDIC) = zeros([size(dataIn2,1) ...
                                             size(dataIn2,2) ...
                                             numel(indexDIC)]);
            end
            if  (max(indexGreen) > 0 && max(indexRed) > 0)
                dataL(:,:,indexRed) = dataL_Red+numNeutropGreen;
                dataL(dataL==numNeutropGreen) = 0;
                dataL(:,:,indexGreen) = dataL_green;
                numNeutrop = numNeutropGreen+numNeutropRed;
            elseif  (max(indexRed) > 0)
                dataL(:,:,indexRed) = dataL_Red;
                numNeutrop = numNeutropRed;            
            elseif  (max(indexGreen) > 0)
                dataL(:,:,indexGreen) = dataL_green;
                numNeutrop = numNeutropGreen;     
            else
                dataL = zeros(handles.rows,handles.cols,handles.levs);
                numNeutrop = 0;     
           end

        if handles.minBlob~=0

            statsData = regionprops(dataL); 
            levsDataL = size(dataL,3);
            %first, find all those that are smaller than 2*handles.minBlob, 
            % merge if they nearly overlap with another
            %
%             currentVolumes                          =([statsData.Area]<(2*handles.minBlob));
%             if any(currentVolumes)
%                 % find the candidate to merge and the overlaps                
%                 for counterneut =1:numNeutrop
%                     if currentVolumes(counterneut)==1
%                         currentCandBB                   = statsData(counterneut).BoundingBox;
%                         if numel(currentCandBB)==4
%                             currentCandBB               = [currentCandBB(1) currentCandBB(2) 1 currentCandBB(3) currentCandBB(4)  0];
%                         end
%                         currentcandidatex               = max(1,floor(currentCandBB(1))):max(1,floor(currentCandBB(1)))+currentCandBB(4);
%                         currentcandidatey               = max(1,floor(currentCandBB(2))):max(1,floor(currentCandBB(2)))+currentCandBB(5);
%                         currentcandidatez               = max(1,floor(currentCandBB(3))):min(levsDataL,max(1,floor(currentCandBB(3)))+currentCandBB(6));
%                         for counterneut2=1:numNeutrop
%                             if counterneut2~=counterneut
%                                 possibleParentBB        =statsData(counterneut2).BoundingBox;
%                                 if numel(possibleParentBB)==4
%                                     possibleParentBB               = [possibleParentBB(1) possibleParentBB(2) 1 possibleParentBB(3) possibleParentBB(4)  0];
%                                 end
% 
%                                 possibleparentx         =max(1,floor(possibleParentBB(1))):max(1,floor(possibleParentBB(1)))+possibleParentBB(4);
%                                 possibleparenty         =max(1,floor(possibleParentBB(2))):max(1,floor(possibleParentBB(2)))+possibleParentBB(5);
%                                 possibleparentz         =max(1,floor(possibleParentBB(3))):min(levsDataL,max(1,floor(possibleParentBB(3)))+possibleParentBB(6));
% 
%                                 overlapInz              = (intersect(possibleparentz,currentcandidatez));
%                                 distanceInX             =min(([min(abs(currentcandidatex(1)-possibleparentx)) min(abs(currentcandidatex(end)-possibleparentx)) ]));
%                                 distanceInY             =min(([min(abs(currentcandidatey(1)-possibleparenty)) min(abs(currentcandidatey(end)-possibleparenty)) ]));
%                                 distanceInZ             =min(abs([possibleParentBB(1)-currentCandBB(1) possibleParentBB(1)-(currentCandBB(1)+currentCandBB(4)) (possibleParentBB(1)+possibleParentBB(4))-currentCandBB(1) (possibleParentBB(1)+possibleParentBB(4))-(currentCandBB(1)+currentCandBB(4))]));
% 
%                                 %if the distance between candidates is <2, AND There is overlap in Z then join
%                                 if ((distanceInX+distanceInY)<3)&&(~isempty(overlapInz))
%                                     %change the class of the currentcandidate
% 
%                                     %bridge between them
%                                     totalRegionX =min(possibleparentx(1),currentcandidatex(1)):max(possibleparentx(end),currentcandidatex(end));
%                                     totalRegionY =min(possibleparenty(1),currentcandidatey(1)):max(possibleparenty(end),currentcandidatey(end));
%                                     for counterLev = currentcandidatez
% 
%                                         %try
%                                         currDataL   = dataL(totalRegionY,totalRegionX,counterLev);
%                                         %catch
%                                         %    qq=1;
%                                         %end
%                                         tempDataL1 = (bwdist(currDataL==counterneut));
%                                         tempDataL2 = (bwdist(currDataL==counterneut2));
%                                         tempDataL3 = ((tempDataL1+tempDataL2)<=3).*(currDataL==0);
%                                         if max(tempDataL3(:))>0
%                                             %try
%                                             dataL(totalRegionY,totalRegionX,counterLev)= dataL(totalRegionY,totalRegionX,counterLev)+(tempDataL3)*counterneut2;
%                                             %catch
%                                             %    qqq=1;
%                                             %end
%                                         end
%                                     end
%                                     dataL(dataL==counterneut) = counterneut2;
% 
%                                 end                               
%                             end
%                         end
%                     end
%                 end    
% 
%             end

            regS = regionprops(dataL,'area');

            %second delete all those that are smaller than handles.minBlob

            if (max(indexRed)>0)&&(max(indexGreen)>0)
                %in case there are two channels, the final labelling has to be separate
                [dataL_green,numNeutropGreen] = ...
                    bwlabeln(ismember(...
                            dataL(:,:,indexGreen),...
                            find([regS.Area]>handles.minBlob(1))...
                            ));
                [dataL_red,numNeutropRed] = ...
                    bwlabeln(ismember(...
                            dataL(:,:,indexRed),...
                            find([regS.Area]>handles.minBlob(2))...
                            ));

                dataL(:,:,indexGreen) = dataL_green;
                dataL(:,:,indexRed) = dataL_red+max(dataL_green(:)).*(dataL_red>0);                    
                statsData = regionprops(dataL);
                numNeutrop = numNeutropGreen+numNeutropRed;

            else
                [dataL,numNeutrop] = ...
                    bwlabeln(ismember(...
                            dataL,...
                            find([regS.Area]>handles.minBlob)...
                            ));
                statsData = regionprops(dataL);
            end
                
        else
            statsData = regionprops(dataL,'area'); 
        end

        save(dataOutName1,'dataL','numNeutrop','statsData');
    end

    dataL = dataOutName;
else
	% --------------------------- dataIn is a Matlab Matrix ---------------------------

    % handles.thresLevel can be a single value or a pair of values (high, low)
    % if it is one, split into two to have separate thresholds
    if numel(handles.thresLevel)==1
        lowT_R = 0.95*(handles.thresLevel);
        highT_R = 1.05*(handles.thresLevel);            
        lowT_G = lowT_R;
        highT_G = highT_R;
    elseif numel(handles.thresLevel)==2
        lowT_R = min(handles.thresLevel);
        highT_R = max(handles.thresLevel);      
        lowT_G = lowT_R;
        highT_G = highT_R;
    elseif numel(handles.thresLevel)==4
        lowT_G = min(handles.thresLevel(1:2));
        highT_G = max(handles.thresLevel(1:2));            
        lowT_R = min(handles.thresLevel(3:4));
        highT_R = max(handles.thresLevel(3:4));                      
    else 
        handles=[]; return;
    end
        
    % handles.ChannelDistribution is used to determine which slices 
    % correspond to Green/Red/DIC
  
     if max(indexGreen)>0
        % Only process if there are GREEN channel slices
         [GreenLowT,numGreenLowT] = bwlabeln(dataIn(:,:,indexGreen)>lowT_G);
         if numGreenLowT>0
             GreenHighT = (dataIn(:,:,indexGreen)>highT_G);
             prodHighLowT = GreenLowT.*GreenHighT;
             GreenLowT_2D = (sum(GreenLowT,3))>0;
             clear rr cc
             [rr,cc] = find(GreenLowT_2D);
             rr = unique(rr);
             cc = unique(cc);
             
             neutropToKeep = unique(prodHighLowT(rr,cc,:));
             
             [dataL_green,numNeutropGreen] = ...
                 bwlabeln(ismember(GreenLowT,neutropToKeep(2:end))); 
         else
             %This line assumes that there is only green channel
             %dataL_green                        = zeros(size(dataIn));
             %Fixed here:
             dataL_green = zeros([size(dataIn,1) size(dataIn,2) numel(indexGreen)]);
             numNeutropGreen = 0;
         end
     else
         numNeutropGreen = 0;
     end

     if max(indexRed)>0
         % Only process if there are RED channel slices
         [RedLowT,numRedLowT] = bwlabeln(dataIn(:,:,indexRed)>lowT_R);
         if numRedLowT>0
             RedHighT = (dataIn(:,:,indexRed)>highT_R);
             prodHighLowT = RedLowT.*RedHighT;
             RedLowT_2D = (sum(RedLowT,3))>0;
             clear rr cc
             [rr,cc] = find(RedLowT_2D);
             rr = unique(rr);
             cc = unique(cc);
             
             neutropToKeep = unique(prodHighLowT(rr,cc,:));
             
             [dataL_Red,numNeutropRed] = ...
                 bwlabeln(ismember(RedLowT,neutropToKeep(2:end)));                   
         else
             %This line assumes that there is no green channel, only red.
             %dataL_Red                          = zeros(size(dataIn2));
             %Fixed here:
             dataL_Red = zeros([size(dataIn,1) size(dataIn,2) numel(indexRed)]);
             numNeutropRed = 0;
         end
     else
         numNeutropRed = 0;
     end

    % now merge BACK into a single matrix if only one case
    if (max(indexDIC) > 0)
        dataL(:,:,indexDIC) = ...
            zeros([size(dataIn,1) size(dataIn,2) numel(indexDIC)]);
    end
    if ( max(indexGreen)>0 && max(indexRed)>0 )
        dataL(:,:,indexRed) = dataL_Red+numNeutropGreen;
        dataL(dataL==numNeutropGreen) = 0;
        dataL(:,:,indexGreen) = dataL_green;
        numNeutrop = numNeutropGreen+numNeutropRed;

    elseif  (max(indexRed)>0)
        dataL = dataL_Red;
        numNeutrop  = numNeutropRed;            
    elseif  (max(indexGreen)>0)
        dataL = dataL_green;
        numNeutrop = numNeutropGreen;     
    else
        dataL = zeros(handles.rows,handles.cols,handles.levs);
        numNeutrop = 0;     
   end

    if handles.minBlob~=0
        % removal of objects smaller than a predefined blob level, can be bypassed with handles.minBlob=0 if desired
            regS = regionprops(dataL,'area');
            % delete all those that are smaller than handles.minBlob

            if (max(indexRed)>0)&&(max(indexGreen)>0)
                %in case there are two channels, the final labelling has to be separate
                [dataL_green,numNeutropGreen] = ...
                    bwlabeln(ismember(dataL(:,:,indexGreen),...
                            find([regS.Area]>handles.minBlob)...
                            ));
                
                [dataL_red,numNeutropRed] = ...
                    bwlabeln(ismember(dataL(:,:,indexRed),...
                    find([regS.Area]>handles.minBlob)...
                    ));

                dataL(:,:,indexGreen) = dataL_green;
                dataL(:,:,indexRed) = dataL_red+max(dataL_green(:)).*(dataL_red>0);                    
                statsData = regionprops(dataL);
                numNeutrop = numNeutropGreen+numNeutropRed;
            else
                [dataL,numNeutrop] = ...
                    bwlabeln(ismember(dataL,...
                            find([regS.Area]>handles.minBlob)...
                            ));
                
                statsData = regionprops(dataL);
            end   
    end
end
