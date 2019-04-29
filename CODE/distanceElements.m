function [miuD,distanceBetPoints,numNeighboursAtDist,compDistance,distToDisappearing,distToAppearing] =distanceElements (dataIn)
%function [miuD,distanceBetPoints,numNeighboursAtDist,compDistance,distToDisappearing,distToAppearing] =distanceElements (dataIn)
%
%--------------------------------------------------------------------------
% distanceElements  obtains distance between a series of elements, it returns 
%     the minimum distance from the element to the neighbours, i.e. the 
%     distance to the closest neighbour
%
%       INPUT
%         dataIn:               An matrix of x,y,z positions or a 2D BW image.
%       OUTPUT
%         miuD:                 Distance to other points as integer (ceil)
%         distanceBetPoints:    Distance to other points
%         numNeighboursAtDist:  Count of neighbours at certain distance,
%                               number is coded.
%         comDistance:          Distance to closest neighbour
%         distToDisappearing:   Distance to closest disappearing neighbour
%         distToAppearing:      Distance to closest appearing neighbour
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

if isa(dataIn,'logical')
    dataIn=bwlabel(dataIn);
    centroidsIn=regionprops(dataIn,'Centroid');
    dataIn=cell2mat(struct2cell(centroidsIn)');
end
if isstruct (dataIn)
    dataIn=cell2mat(struct2cell(dataIn)');
end


[rows,cols,levs]= size(dataIn);

if rows==1
    miuD                        =inf;
    distanceBetPoints           =[];
    numNeighboursAtDist         =0;
    compDistance                =[];
        distToDisappearing      =inf;
    distToAppearing             =inf;

else
    distBetPoints               =[];
    numNeighboursAtDist         =[];
    compDistance                =[];
    distToDisappearing          =[];
    distToAppearing             =[];
    for counterR=1:rows
        indexRows=[1:counterR-1  counterR+1:rows];
        %----- finds distances from every element to all the rest
        diffDistances           = (dataIn(indexRows,1:3)-repmat(dataIn(counterR,1:3),[ rows-1 1 ]));
        absDistances            = sqrt( sum((diffDistances).^2,2));
        [minAbsD,indexAbsD]     = min(absDistances);
        %----- counts how many neighbours at certain distance brackets and stores every unit of 10
        numNeighbours           =   (    sum((absDistances>0)&(absDistances<=10))...
                                    +10* sum((absDistances>10)&(absDistances<=20))...
                                    +1e2*sum((absDistances>20)&(absDistances<=30))...
                                    +1e3*sum((absDistances>30)&(absDistances<=40))...
                                    +1e4*sum((absDistances>40)&(absDistances<=50))...
                                    +1e5*sum((absDistances>50)&(absDistances<=100))...
                                    +1e6*sum((absDistances>100))) ;
        %----- number of neighbours                    
        numNeighboursAtDist     = [numNeighboursAtDist;numNeighbours];
        %----- distance to closest neighbour in absolute distance
        distBetPoints           = [distBetPoints;minAbsD];
        %----- distance to closest neighbour in x,y
        compDistance            = [compDistance;diffDistances(indexAbsD,1)+diffDistances(indexAbsD,2)*i];
        %----- distance to closest DISAPPEARING NEIGHBOUR only possible if data of Parent/Child is
        %      provided [X Y Z Parent (7) Child (8) ], it is provided in absolute terms
        if (cols ==5)
            distToDisappearingTemp                  = absDistances;
            distToAppearingTemp                     = absDistances;
            distToDisappearingTemp(dataIn(indexRows,5)~=0)  = inf;
            distToAppearingTemp(dataIn(indexRows,4)~=0)     = inf;
            distToDisappearing  = [distToDisappearing;min(distToDisappearingTemp)];
            distToAppearing     = [distToAppearing;   min(distToAppearingTemp)];
            
        else
            distToDisappearing  = [distToDisappearing;inf];
            distToAppearing     = [distToAppearing;inf];
            
        end
        
    end
    
    distanceBetPoints=(distBetPoints);
    miuD=ceil(distanceBetPoints);
end
