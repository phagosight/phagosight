function f=drawsquare_r(xx1,yy1,xx2,yy2,lineW,dlt,maxz,axestag,viewax)
%function f=drawsquare_r(xx1,yy1,xx2,yy2,lineW,dlt,maxz,axestag,viewax)
%
%--------------------------------------------------------------------------
% drawsquare_r  draws a square of white and black lines over the current axes
%
%       INPUT   
%         xx1-xx2:      dimensions in columns, yy1-yy2 dimensions in rows
%         lineW:        width of the primary line       (white)
%         dlt:          separation with secondary line  (black)
%         maxz:         verifies that the level of display is adequate for the data
%         axestag:      tags for the object (they can later be manipulated)
%         viewax:       axis on which the image is going to be printed        
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


%--------------------------------------------------------------------------
%careful ****** x - column   y -row!!!!!!!!!!!!!
%------ no input data is received, error -------------------------
if nargin < 4
    help drawsquare_r;  return;
end;


%------ only position is received ------
if (nargin<=5)
    %----- get the axis on which the image is going to be printed
    viewax=gca;
    %----- get the z axis limits to determine the value to set the height
    zlimits=get(viewax,'Zlim');
    %------ verify size and set in correct order
    x1=max(min(xx1,xx2),1);
    y1=max(min(yy1,yy2),1);
    x2=max(xx1,xx2);
    y2=max(yy1,yy2);
    maxz=zlimits(1)+1;
    axestag='t_alin';
    if nargin==4
        lineW=4;
    end;
    dlt=lineW/2;
end
%set(figure(1),'windowbuttonmotionfcn',['f=round(get(gca,''currentpoint''));','drawlines(i,f)'])
if ~exist('viewax','var'); viewax=gca; end
if ~exist('x1','var') 
    x1=max(min(xx1,xx2),1); 
    y1=max(min(yy1,yy2),1);
    x2=max(xx1,xx2);
    y2=max(yy1,yy2);
end
if ~exist('maxz','var');        zlimits=get(viewax,'Zlim'); maxz=zlimits(1)+1;  end;
if ~exist('axestag','var');     axestag='t_alin';                               end
if ~exist('lineW','var');       lineW=4;                                        end
if ~exist('dlt','var');         dlt=lineW/2;                                    end



line([x1 x1],[y1 y2],[maxz maxz],'color',[1 1 1],'LineWidth',lineW,'Tag',strcat(axestag,'1'),'Parent',viewax)
line([x1 x2],[y1 y1],[maxz maxz],'color',[1 1 1],'LineWidth',lineW,'Tag',strcat(axestag,'2'),'Parent',viewax);
line([x2 x2],[y1 y2],[maxz maxz],'color',[1 1 1],'LineWidth',lineW,'Tag',strcat(axestag,'3'),'Parent',viewax);
line([x1 x2],[y2 y2],[maxz maxz],'color',[1 1 1],'LineWidth',lineW,'Tag',strcat(axestag,'4'),'Parent',viewax);
line([x1+dlt x1+dlt],[y1+dlt y2-dlt],[maxz maxz],'color',[0.02 0 0],'LineWidth',lineW/4,'Tag',strcat(axestag,'5'),'Parent',viewax);
line([x1+dlt x2-dlt],[y1+dlt y1+dlt],[maxz maxz],'color',[0.02 0 0],'LineWidth',lineW/4,'Tag',strcat(axestag,'6'),'Parent',viewax);
line([x2-dlt x2-dlt],[y1+dlt y2-dlt],[maxz maxz],'color',[0.02 0 0],'LineWidth',lineW/4,'Tag',strcat(axestag,'7'),'Parent',viewax);
line([x1+dlt x2-dlt],[y2-dlt y2-dlt],[maxz maxz],'color',[0.02 0 0],'LineWidth',lineW/4,'Tag',strcat(axestag,'8'),'Parent',viewax);


