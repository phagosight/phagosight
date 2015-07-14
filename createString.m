function [str] = createString(pxi,pyi,pzi,pti,pxf,pyf,pzf,ptf,vx,vy,vz,vt,vi,mIn,mIn2,mIn4,mIn3,mIn5)
%function [str] = createString(pxi,pyi,pzi,pti,pxf,pyf,pzf,ptf,vx,vy,vz,vt,vi,mIn,mIn2,mIn4,mIn3,mIn5)
%
%--------------------------------------------------------------------------
% createString  creates a description of the track and vertex information 
%     provided in the arguments.
%
%       INPUT
%           pxi:	track's initial point, x-coordinate.
%           pxf:    track's final point, x-coordinate.
%           pyi:    track's initial point, y-coordinate.
%           pyf:    track's final point, y-coordinate.
%           pzi:    track's initial point, z-coordinate.
%           pzf:    track's final point, z-coordinate.
%           pti:    track's initial point, time point (frame).
%           ptf:    track's final point, time point (frame).
%           vx:     selected vertex, x-coordinate.
%           vy:     selected vertex, y-coordinate.
%           vz:     selected vertex, z-coordinate.
%           vt:     selected vertex, time point (frame).
%           vi:     selected vertex id.
%           mIn:    meander Ratio
%           mIn2:   tortuosity
%           mIn4:   length of Track
%           mIn3:   distance
%           mIn5:   average distance per hop
%
%       OUTPUT  
%           str:    test description.
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
	str = sprintf([' Position  (initial-final):\n  X: \t %d \t %d \n'...
		'  Y: \t %d \t %d\n  Z: \t %d \t %d\n  T: \t %d \t %d  \n\n'...
		'  Vertex: (point selected)\n  X: \t %10.3f\n \tY:   %10.3f\n  '...
		'Z: \t %10.3f\n  T: \t %d  \n\n Track: \t %d\n \n Meander Ratio:'...
		'    \t %10.3f \n \n Tortuosity:           \t%10.3f \n \n '...
		'Length  Track:    \t %10.3f \n \n Dist  (init,fin):     \t %10.3f \n \n'...
        'Av Dist per Frame: \t %10.3f' ],...
		pxi,pxf,pyi,pyf,pzi,pzf,pti,ptf,vx,vy,vz,vt,vi,mIn,mIn2,mIn4,mIn3,mIn5);
    userhandles = get(gcf,'userdata');
    userhandles.trackValue = vi;

    set(gcf,'userdata',userhandles);
