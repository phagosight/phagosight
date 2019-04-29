function bor = zerocross(imag)
%function bor = zerocross(imag)
%
%--------------------------------------------------------------------------
% zerocross  get zero crossings of a certain input can be 2D or 1D
%
%       INPUT
%         imag:     data (image or signal) from which the zero crossing
%                   points are required 
%
%       OUTPUT
%         bor:      matrix with 1 where there is a zero cross  and 0 elsewhere  
%                   same dimensions as the input
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

[lins,cols,levels] = size(imag);
delta = 0.00001;

if ~isa(imag,'logical')
    imag=sign(imag);
end
%------ revise the case ------
%------    1 1D use plot for either line or column data
%------    2 1D but not line or column, stored in various z
%------    3 2D use surfdat_r.m function
%------    4 3D use just the base of cube

if ((cols==1|lins==1)&levels==1)  %- 1D over main plane
    if cols>= lins
        yy = [0 imag(1:cols-1)];
    else
        yy = [0 imag(1:lins-1)']';
    end;
    bor = abs((sign(imag+delta)-sign(yy+delta)));
elseif (cols==1&lins==1&levels~= 1)%- 1D over other plane
    imag = permute(imag,[2 3 1]);
    yy = [0 imag(1:cols-1)];
    bor = (sign(imag+delta)-sign(yy+delta));
elseif (lins~= 1&cols~= 1&levels==1) %- 2D over main plane
    %------ only 1 degree neighbourhood considered---------
    diffVer         = diff(imag,1,1);zerCols = zeros(1,cols);
    diffHor         = diff(imag,1,2);zerRows = zeros(lins,1);
    qq1             = [zerCols;(diffVer)>0];
    qq2             = [(diffVer)<0;zerCols];
    qq3             = [ (diffHor)<0 zerRows ];
    qq4             = [ zerRows (diffHor)>0 ];
    bor             = qq1|qq2|qq3|qq4;
elseif(lins~= 1&cols~= 1&levels~= 1) %- 3D
    yy5             = [zeros(1,cols,levels);  imag(1:lins-1,1:cols,:)];     %|d
    yy6             = [imag(2:lins,1:cols,:);zeros(1,cols,levels)];                  %u|
    yy7             = [ zeros(lins,1,levels) imag(1:lins,1:cols-1,:)];               %-r
    yy8             = [ imag(1:lins,2:cols,:) zeros(lins,1,levels)];                 %l-
    bor5            = fix(delta+(sign(imag)-sign(yy5))/2);
    bor6            = fix(delta+(sign(imag)-sign(yy6))/2);
    bor7            = fix(delta+(sign(imag)-sign(yy7))/2);
    bor8            = fix(delta+(sign(imag)-sign(yy8))/2);
    bor             = sign(bor5+bor6+bor7+bor8);
end;
