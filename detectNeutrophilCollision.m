function collisionNetwork = detectNeutrophilCollision(handles)
%function collisionNetwork = detectNeutrophilCollision(handles)
%
%--------------------------------------------------------------------------
% detectNeutrophilCollision  scans the tracks through time to detect any two
%     neutrophils that may have collided and merged into a single blob
%
%       INPUT  
%         handles:          handles struct
%       OUTPUT  
%         collisionNetwork: handles.finalNetwork updated containing
%             information about collisions
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


%% Define some used values
numMetric                               = 10;
numMetric2                              = 30;
numMetric3                              = 31;
%%

sortedVols                              = (sort(handles.nodeNetwork(:,numMetric)));
numObjects                              = size(handles.nodeNetwork,1);
indexNumObjects                         = (1:numObjects)/numObjects;


minObjVolume                            = sortedVols(find(indexNumObjects>0.05,1,'first'));
medObjVolume                            = sortedVols(find(indexNumObjects>0.5,1,'first'));
maxObjVolume                            = sortedVols(find(indexNumObjects>0.90,1,'first'));

% the distance to the disappearing object is determined by the radius of the medium object
radMedObj                               = (0.75*medObjVolume/pi)^(1/3);
outlierDis                              = floor(radMedObj*2);

numTrackNet                               = size( handles.finalNetwork,2);
collisionNetwork                        = handles.finalNetwork;
for numTrack = 1:numTrackNet
    % Proceed track by track
    currTrack                               = handles.finalNetwork(:,numTrack);
    [maxR]                                  = find(currTrack);

    currTrackEff                            = currTrack(maxR); %#ok<FNDSB>
    finNode                                 = size(currTrackEff,1);

    trackVol                                = handles.nodeNetwork(currTrackEff,numMetric);
    trackDisDis                             = handles.nodeNetwork(currTrackEff,numMetric2);
    trackDisApp                             = handles.nodeNetwork(currTrackEff,numMetric3);

    neigh0_10   = rem(handles.nodeNetwork(currTrackEff,28),10);
    neigh10_20  = rem(floor(handles.nodeNetwork(currTrackEff,28)/10), 10);
    neigh20_30  = rem(floor(handles.nodeNetwork(currTrackEff,28)/100),10);

    trackNum                                = neigh0_10+neigh10_20+neigh20_30;

    diffVol                                 = diff(trackVol);


    %
    % here the rules for a  * M E R G E R * are defined (H -Hard rule,  needs to happen S-soft rule, may happen)
    % REDEFINED after distance to disappearing/appearing was determined:
    % 1 H - There should be an increase in volume, when 2 of equal size merge or a small with a big is easy,
    %       when one big eats a small one, change is relatively small. Take a change > 5% of sizes of objects
    % 2 H - there should be at least one neighbour which disappeared in the past frame, and this should be closer than 20 pix
    %       and whose volume is similar to the increase


    RuleMerge1                              = [0;diffVol     >  minObjVolume];
    RuleMerge2                              = [0;trackDisDis(1:end-1) < outlierDis];

    RuleMerge3                              = trackVol    >  maxObjVolume;
    RuleMerge4                              = [0;diffVol  >  medObjVolume];
    RuleMerge5                              = [0;trackNum(1:end-1)]    >  0;

    
    DefMerger1                               = (RuleMerge1&RuleMerge2);
    DefMerger2                              = (RuleMerge3&RuleMerge4&RuleMerge5);
    DefMergerT                              = (DefMerger1&DefMerger2);
    %discard any changes that happen in the first or last frame
    DefMergerT([1:2 end-1:end])             =0;
    DefMergerC                              = cumsum(DefMergerT);

    
    
    % here the rules for a  * SPLIT * are defined (H -Hard rule,  needs to happen S-soft rule, may happen)
    % REDEFINED after distance to disappearing/appearing was determined:
    % 1 H - There should be an decrease in volume, when 2 of equal size merge or a small with a big is easy,
    %       when one big eats a small one, change is relatively small. Take a change > 5% of sizes of objects
    % 2 H - there should be at least one neighbour which appeared in the THAT frame, and this should be closer than 20 pix
    %       and whose volume is similar to the increase
    % Split can be confused by neighbouring objects that do not reflect rules 1 & 2 so use an alternative
    % Harder Rule:
    % 3 S - The original Volume must be Big for a Start, larger than 85% of objects
    % 4 S - There must be a big decrease in volume, not just of the min object (5%) but of 50%
    % 5 S - There must be close neighbours (that introduced the confusion for a start)
    %RuleSplit1                              = diffVol     <-  minObjVolume;
    RuleSplit1                              = [0;diffVol     <-  minObjVolume];
    RuleSplit2                              = trackDisApp < outlierDis;

    RuleSplit3                              = trackVol    >  maxObjVolume;
    RuleSplit4                              = [0;diffVol     <-  medObjVolume];
    RuleSplit5                              = trackNum    >  0;


    DefSplit1       = (RuleSplit1&RuleSplit2);

    DefSplit2       = (RuleSplit3&RuleSplit4&RuleSplit5);
    DefSplitT       = (DefSplit1&DefSplit2);
    
    %discard any changes that happen in the first or last frame
    DefSplitT([1:2 end-1:end])              = 0;

    DefSplitC       = cumsum(DefSplitT);
    CumulObj        = DefMergerC-DefSplitC;
    %% To avoid having a fluke split of 1 or 2 frames, discard all those split before mergers (negative) <2 frames
    [negAreasSplit]             = find(CumulObj==-1);
    [StartSplits]               = find((DefSplitT));

    negAreasStartSplit          = StartSplits(ismember(StartSplits,negAreasSplit));
    for counterSplits=1:size(negAreasStartSplit)
        startCount = max(1,negAreasStartSplit(counterSplits)-2);
        endCount   =  min(finNode,negAreasStartSplit(counterSplits)+3);
        
        if sum(CumulObj(startCount:endCount)==-1)<4
            CumulObj(startCount:endCount)=0;
        end
    end
    %%
    
    CumulObj        =1+(CumulObj)- min(CumulObj);
    
    collisionNetwork(1:finNode,numTrack) =CumulObj;

end


