function selectNeutrophilsM(arg,handles,woundRegion)
%function selectNeutrophilsM(arg,handles,woundRegion)
%--------------------------------------------------------------------------
% selectNeutrophilsM  is a GUI based function that selects tracks and displays 
%     some of their characteristics it is necessary to pass a figure where
%     the tracks have been displayed through the use of plotTracks. This file 
%     is based on select3dtool, downloaded from matlab central
%
%       INPUT
%         arg:          the current axis where the tracks have been displayed
%         handles:      handles struct with details of the tracks
%         woundRegion:  wound region matrix in case tracks are recalculated
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
% This m-file is part of the PhagoSight package used to analyse
% fluorescent phagocytes as observed through confocal or
% multiphoton microscopes. For a comprehensive 
% user manual, please visit:
%
%           <http://www.phagosight.org.uk>
%
% Please feel welcome to use, adapt or modify the files. If you can improve
% the performance of any other algorithm please contact us so that we can
% update the package accordingly.
%
%--------------------------------------------------------------------------
%
% The authors shall not be liable for any errors or responsibility for the 
% accuracy, completeness, or usefulness of any information, or
% method in the content, or for any actions taken in reliance thereon.
%
%--------------------------------------------------------------------------
    
    if (nargin<1)||(isempty(arg))
        arg = gcf;
    end
    
    if ~ishandle(arg)
        feval(arg);
        return;
    end
    %%

    if ~isfield(handles.distanceNetwork,'meanderRatio')
        if ~exist('woundRegion','var')
            help selectNeutrophilsM
            return;
        else
            handles=effectiveDistance(handles,woundRegion);
            handles=effectiveTracks(handles,woundRegion);
        end
    end


    if ~isfield(handles.distanceNetwork,'tortuosity')
            % Tortuosity and meanderRatio are inverse of each other
         handles.distanceNetwork.tortuosity = ...
             1./handles.distanceNetwork.meanderRatio;
    end

    %% initialize gui %%
    fig = arg;
    figure(fig);

    uistate = uiclearmode(fig);
    [tool, htext] = createUI;
    
    v=ver('matlab');
    if str2num(v.Release(3:6)) < 2014
        hmarker1 = line('marker','o','markersize',10,'markerfacecolor',...
            'k','erasemode','xor','visible','off');
        hmarker2 = line('marker','o','markersize',10,'markerfacecolor',...
            'r','erasemode','xor','visible','off');
    else
        hmarker1 = line('marker','o','markersize',10,'markerfacecolor',...
            'k','visible','off');
        hmarker2 = line('marker','o','markersize',10,'markerfacecolor',...
            'r','visible','off');
    end

    state.uistate = uistate;
    state.text = htext;
    state.tool = tool;
    state.fig = fig;
    state.marker1 = hmarker1;
    state.marker2 = hmarker2;
    setappdata(fig,'selectNeutrophilsM',state);
    setappdata(state.tool,'select3dhost',fig);

    set(fig,'windowbuttondownfcn','selectNeutrophilsM(''click'')');


    if exist('handles','var')
        userhandles.handles         = handles;
    end

    if exist('woundRegion','var')
        userhandles.woundRegion     = woundRegion;
    end
    userhandles.trackValue          = 0;
    
    set(fig,'userdata',userhandles);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function off
    
    state = getappdata(gcbf,'selectNeutrophilsM');
    
    if ~isempty(state)
        delete(state.tool);
    end
    
    fig = getappdata(gcbf,'select3dhost');
    
    if ~isempty(fig) & ishandle(fig)
        state = getappdata(fig,'selectNeutrophilsM');     
        uirestore(state.uistate);
        delete(state.marker1);
        delete(state.marker2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function click
    
    [p v vi] = select3d;
    state = getappdata(gcbf,'selectNeutrophilsM');

    if ~ishandle(state.text)
        state.text = createUI;
    end

    if ~ishandle(state.marker1)
        state.marker1 = [];
    end

    if ~ishandle(state.marker2)
        state.marker2 = [];
    end
    
    setappdata(state.fig,'selectNeutrophilsM',state);
    
    if isempty(v)
        v = [nan nan nan nan];
        vi = nan;
        set(state.marker2,'visible','off');
    else
        userhandles         = get(gcf,'userdata');
        %%
        handles = userhandles.handles;
        currFrame = handles.nodeNetwork(handles.nodeNetwork(:,5)==v(3),:);
        distToSelected = currFrame(:,[2 1])-repmat(v(1:2),[size(currFrame,1) 1]);
        [minD,indMinD] = min(sqrt(sum(distToSelected.^2,2)));
        trackSelected = currFrame(indMinD,13);
        if trackSelected==0
            currTrack = handles.nodeNetwork(currFrame(indMinD,6),:);
            mIn = 0;
            mIn2 = 0;
            mIn4 = 0;
            mIn3 = 0;
            mIn5 = 0;
        else
            currTrack = handles.nodeNetwork(...
                handles.nodeNetwork(:,13)==trackSelected,:);
            mIn = handles.distanceNetwork.meanderRatio(trackSelected);
            mIn2 = handles.distanceNetwork.totPerTrack(trackSelected);
            mIn4 = handles.distanceNetwork.tortuosity(trackSelected);
            mIn3 = mIn2/mIn4;
            mIn5 = handles.distanceNetwork.avPerTrack(trackSelected);
            
        end
        p = round( [ currTrack(1,[2 1 3 5]) currTrack(end,[2 1 3 5])]);

        v2 = [v(1) v(2) currTrack(currTrack(:,5)==v(3),3)  v(3)];
        vi = trackSelected;



        %%    
        set(state.marker2,'visible','on','xdata',v(1),'ydata',v(2),'zdata',v(3));
    end

    if isempty(p)
        p = [nan nan nan nan];
        set(state.marker1,'visible','off');
    else
        set(state.marker1,'visible','on','xdata',p(1),'ydata',p(2),'zdata',p(4));
    end
    % Update tool and markers
    if numel(p)==8
        set(state.text,'string',createString(...
            p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),...
            v2(1),v2(2),v2(3),v2(4),vi,...
            mIn,mIn4,mIn2,mIn3,mIn5));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig, h] = createUI
    
    pos = [100 500 250 350];

    % Create selection tool %
    fig = figure('handlevisibility','off','menubar','none','resize','off',...
                 'numbertitle','off','name','Select 3-D Tool','position',...
                 pos,'deletefcn','selectNeutrophilsM(''off'')');

    h = uicontrol('style','text','parent',fig,'string',...
                  createString(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),...
                  'units','norm','position',[0 0 1 1],'horizontalalignment','left');
    
    deleteTrackCB = ['userhandles = get(gcf,''userdata'');',...
                     'if userhandles.trackValue~=0;',...
                     'if (isfield(userhandles,''woundRegion''));',...
        'handles=removeOneTrack(userhandles.handles,userhandles.trackValue,userhandles.woundRegion);',...
        'else;handles=removeOneTrack(userhandles.handles,userhandles.trackValue);',...
        'end;end;',...
        'plotTracks(handles.nodeNetwork(:,[1 2 5]), handles.finalNetwork(:,: ),2,1,1);axis([1 handles.cols 1 handles.rows 1 handles.numFrames]);',...
        'userhandles.handles= handles; userhandles.trackValue          = 0;',...
        'set(gcf,''userdata'',userhandles);view(0,0);'];

