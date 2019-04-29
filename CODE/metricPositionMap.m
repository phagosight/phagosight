function [metricMap,metricMap2,B_mean] = metricPositionMap(tempHandles)
%function [metricMap,metricMap2,B_mean] = metricPositionMap(tempHandles)
%
%--------------------------------------------------------------------------
% metricPositionMap  a function that receives the handles with x,y and a 
%     metric and returns a square map with the average value of the metric 
%     per position
%
%       INPUT
%         tempHandles:          handles struct
%
%       OUTPUT
%         metricMap:            metric map        
%         metricMap2:           metric map
%         B_mean:               
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

%% find maximum values of coordinates to places in a 256x 512x 1024 grid  Reduce the grid by 2 to speed up process
maxValueRC                              = max( max(tempHandles(:,1:2)));
sizeGrid                                = 2^ceil(log2(maxValueRC))/2;
[numObjects,numMetrics]                 = size(tempHandles);
%% create the indices that will be used to position pixel locations

tempHandles2(:,4)                        = sub2ind([sizeGrid sizeGrid],round(tempHandles(:,1)/2),round(tempHandles(:,2)/2));
tempHandles2(:,5)                        = 1:numObjects;
indexPositions                          = sort(unique(tempHandles2(:,4)));
numUniquePositions                      = size(indexPositions,1);

B_mean(sizeGrid,sizeGrid,numMetrics-2)  = 0;
B_mean2(2*sizeGrid,2*sizeGrid,numMetrics-2)  = 0;
metricMap(2*sizeGrid,2*sizeGrid,numMetrics-2)  = 0;

B_max(sizeGrid,sizeGrid,numMetrics-2)  = 0;
B_max2(2*sizeGrid,2*sizeGrid,numMetrics-2)  = 0;
metricMap2(2*sizeGrid,2*sizeGrid,numMetrics-2)  = 0;


%% loop over the array to position values
BBB2(1:sizeGrid*sizeGrid,numMetrics-1)=1:sizeGrid*sizeGrid;
BBB3(1:sizeGrid*sizeGrid,numMetrics-1)=1:sizeGrid*sizeGrid;
for k=1:numUniquePositions

    indexN                              = tempHandles2(tempHandles2(:,4)==indexPositions(k),5);

    if ~isempty(indexN)
        for counterMetric = 3:numMetrics
            tempData                    = tempHandles(indexN,counterMetric);

            BBB2(indexPositions(k),counterMetric-2)   = mean(tempData);
            BBB3(indexPositions(k),counterMetric-2)   = sign(mean(tempData))*max(abs(tempData));
        end
    end
end

%%

for counterMetric = 3:numMetrics
    B_mean(:,:,counterMetric-2)  = reshape(BBB2(:,counterMetric-2),[sizeGrid sizeGrid]);
    B_max(:,:,counterMetric-2)   = reshape(BBB3(:,counterMetric-2),[sizeGrid sizeGrid]);
end
B_mean(isnan(B_mean)) =0;
B_max(isnan(B_max)) =0;

filtG = gaussF(5,5,1);

for counterMetric = 3:numMetrics

    B_mean2(:,:,counterMetric-2) = expandu(B_mean(:,:,counterMetric-2));
    metricMap(:,:,counterMetric-2) = imfilter(B_mean2(:,:,counterMetric-2),filtG) ;
    B_max2(:,:,counterMetric-2) = expandu(B_max(:,:,counterMetric-2));
    metricMap2(:,:,counterMetric-2) = imfilter(B_max2(:,:,counterMetric-2),filtG) ;

end


