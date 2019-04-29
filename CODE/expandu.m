function [expData]=expandu(data,numExpansions)
% function [expData]=expandu(data,numExpansions)
%
%--------------------------------------------------------------------------
% reduceu  function to expand data in uniform levels receives an image (or 
%     stack of images) and expands it in quadtree, expansion is ALWAYS in 
%     factor of 2.
%
%       INPUT
%         data:             data to be reduced. 
%         numExpansions:  	number of reductions (1 : 1000x1000 <- 500X500, 
%                           2: 1000x1000 <- 250x250, etc)
%
%       OUTPUT
%         expData:          the data that has been expanded.
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


%------ no input data is received, error -------------------------
if nargin<1
    help expandu;
    expData=[];
    return;
end;
if isempty(data)
    expData=[];
    return;
end;

if ~exist('numExpansions','var')  numExpansions=1; end

if numExpansions>0
    if numExpansions>1
        data=expandu(data,numExpansions-1);    
    end
    [rows,cols,levels]=size(data);
    

        if levels==1 expData=zeros(rows*2,cols*2); else  expData=zeros(rows*2,cols*2,levels*2);   end
        expData(1:2:end,1:2:end,1:levels)=data;
        expData(1:2:end,2:2:end,1:levels)=data;
        expData(2:2:end,1:2:end,1:levels)=data;
        
        expData(2:2:end,2:2:end,1:levels)=data;
        if levels>1 
            expData(:,:,1:2:end)=expData(:,:,1:levels); 
            expData(:,:,2:2:end)=expData(:,:,1:2:end); 
        end

        
end
