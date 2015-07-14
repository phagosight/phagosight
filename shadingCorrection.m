function [dataOut,errSurface,errSurfaceMin,errSurfaceMax] = shadingCorrection(dataIn,numScales)
%function [dataOut,errSurface,errSurfaceMin,errSurfaceMax] = shadingCorrection(dataIn)
%
%-------- this function corrects the shading from images -----------------------------
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             Research Fellow  Sheffield University                   ----------
%------             http://carlos-reyes.staff.shef.ac.uk                    ----------
%------  5 February 2009                                   ---------------------------
%-------------------------------------------------------------------------------------
% A function to correct the shading (inhomogeneous background intensities) of images.
% The function estimates the background from an envelope of the data and then substracts
% it from the data to return an image without shading. The only assumption is that the
% objects are either darker or brigther than the background (but not both)

% input data:       dataIn          : an image with cells vessels or any other kind of objects
% output data:      dataOut         : image with a uniform background
%                   errSurface      : the shading surface
%                   errSurfaceMin   : the lower envelope
%                   errSurfaceMax   : the higher envelope

% Earlier version of the function were colourBackgroundequalise*.m

%% Parse input data
%------ no input data is received, error -------------------------
if nargin <1;     help shadingCorrection; dataOut=[]; return;  end;


%-------- regular size check and determination of first zoom region
[rows,cols,levs]=size(dataIn);
if ~isa(dataIn,'double'); 
    dataIn=double(dataIn); 
    reconvertTo255 =1;
else
    reconvertTo255 =0;
end

%% Define parameters
if ~exist('numScales','var')
    numScales                           = 55;       % maximum distance of analysis
end
stopCriterion                       = 0.05;     % stop criterion difference between steps
% 
% % initialise parameters
% y_scaleMax(rows,cols,numScales)     = 0;
% y_scaleMin(rows,cols,numScales)     = 0;
% dataOut(rows,cols,levs)             = 0;        %#ok<NASGU>
% sizeFR_S                            = 3;        % size of a Gaussian filter to remove noise
% sizeFC_S                            = sizeFR_S;
% filtG_small                         = gaussF(sizeFR_S,sizeFC_S,1); %define the filter


if levs>1
    % initialise parameters
    %y_scaleMax(rows,cols,numScales)     = 0;
    %y_scaleMin(rows,cols,numScales)     = 0;
    dataOut(rows,cols,levs)             = 0;        %#ok<NASGU>
    %sizeFR_S                            = 3;        % size of a Gaussian filter to remove noise
    %sizeFC_S                            = sizeFR_S;
    %filtG_small                         = gaussF(sizeFR_S,sizeFC_S,1); %define the filter

    %in case the data is 3D (a colour image or a volume) subsample to reduce complexity and call same
    %function for each individual channel (level)
    subSampl=1;
    if nargout==4
        for currLev=1:levs
            %[dataOut_S(:,:,currLev),errSurface_S(:,:,currLev)]=shadingCorrection(reduceu(dataIn(:,:,currLev),subSampl/2)); %#ok<AGROW>
            [dataOut_S(:,:,currLev),errSurface_S(:,:,currLev),errSurfaceMin_S(:,:,currLev),errSurfaceMax_S(:,:,currLev)]=shadingCorrection(dataIn(1:subSampl:end,1:subSampl:end,currLev)); %#ok<AGROW>
            %errSurface(:,:,currLev) = expandu(errSurface_S(:,:,currLev),subSampl/2);         %#ok<AGROW>
            %errSurfaceMax(:,:,currLev) = expandu(errSurfaceMax_S(:,:,currLev),subSampl/2);         %#ok<AGROW>
            %errSurfaceMin(:,:,currLev) = expandu(errSurfaceMin_S(:,:,currLev),subSampl/2);         %#ok<AGROW>
        end
        dataOut=dataOut_S;
        errSurface=errSurface_S;
        errSurfaceMax=errSurfaceMax_S;
        errSurfaceMin=errSurfaceMin_S;
        %dataOut                     = dataIn - errSurface(1:rows,1:cols,:);
        errSurfaceMax               = errSurfaceMax(1:rows,1:cols,:);
        errSurfaceMin               = errSurfaceMin(1:rows,1:cols,:);
        dataOut(dataOut>255)        = 255;
        dataOut(dataOut<0)          = 0;
        dataOut                     = (dataOut-min(dataOut(:)));
        dataOut                     = 255*(dataOut/max(dataOut(:)));

    else
        for currLev=1:levs
            %[dataOut_S(:,:,currLev),errSurface_S(:,:,currLev)]=shadingCorrection(reduceu(dataIn(:,:,currLev),subSampl/2)); %#ok<AGROW>
            [dataOut_S(:,:,currLev),errSurface_S(:,:,currLev),errSurfaceMin_S(:,:,currLev),errSurfaceMax_S(:,:,currLev)]=shadingCorrection(dataIn(1:subSampl:end,1:subSampl:end,currLev)); %#ok<AGROW>
            %errSurface(:,:,currLev) = expandu(errSurface_S(:,:,currLev),subSampl/2);         %#ok<AGROW>
        end
        %dataOut = dataIn - errSurface(1:rows,1:cols,:);
        dataOut=dataOut_S;
        errSurface=errSurface_S;
        dataOut(dataOut>255)=255;
        dataOut(dataOut<0)=0;
        %dataOut                     = (dataOut-min(dataOut(:)));
        %dataOut                     = 255*(dataOut/max(dataOut(:)));
    end
else

    % initialise parameters
    % y_scale*** will keep the intermediate steps, no need to keep (r x c x numScales) as even numbers
    % are not used
    y_scaleMax(rows,cols,ceil(numScales/2)) = 0;
    y_scaleMin(rows,cols,ceil(numScales/2)) = 0;
    dataOut(rows,cols,levs)                 = 0;        %#ok<NASGU>
    sizeFR_S                                = 3;        % size of a Gaussian filter to remove noise
    sizeFC_S                                = sizeFR_S;
    filtG_small                             = gaussF(sizeFR_S,sizeFC_S,1); %define the filter

    %---- Low pass filter to reduce the effect of noise
    y1                                      = conv2(padData(dataIn,ceil(sizeFR_S/2)),filtG_small);
    y_LPF                                   = y1(sizeFR_S+1:end-sizeFR_S,sizeFC_S+1:end-sizeFC_S);

    %clear y_scaleM*;

    %---- adjust a surface to the maxima/minima of the data

    for cStep=1:2:numScales
        %disp(cStep)
        %each scale will find the average of the opposite neighbours of
        %a pixel at different degrees of separation
        cStep2                              = ceil(cStep/2);
        y_scaleMax(:,:,cStep2)              = y_LPF;
        y_scaleMin(:,:,cStep2)              = y_LPF;

        %diagonal neighbours of an 8-connectivity
        tempNW                              = y_LPF (1:rows-2*cStep,1:cols-2*cStep);
        tempSW                              = y_LPF (1+2*cStep:rows,1:cols-2*cStep);
        tempSE                              = y_LPF (1+2*cStep:rows,1+2*cStep:cols);
        tempNE                              = y_LPF (1:rows-2*cStep,1+2*cStep:cols);
        %immediate neighbours of a 4-connectivity
        tempN                               = y_LPF (1:rows-2*cStep,1+cStep:cols-cStep);
        tempW                               = y_LPF (1+cStep:rows-cStep,1:cols-2*cStep);
        tempS                               = y_LPF (1+2*cStep:rows,1+cStep:cols-cStep);
        tempE                               = y_LPF (1+cStep:rows-cStep,1+2*cStep:cols);
        clear tempAv*;

        %find averages of opposites and store vertically in a 3D matrix
        tempAv(:,:,1)                       = (tempNE+ tempSW)/2;
        tempAv(:,:,2)                       = (tempNW+ tempSE)/2;
        tempAv(:,:,3)                       = (tempN + tempS) /2;
        tempAv(:,:,4)                       = (tempW + tempE) /2;
        %find minimum and maximum
        tempAvMax                           = max (tempAv,[],3);
        tempAvMin                           = min (tempAv,[],3);
        %at each scale compare the averages with the actual value, keep the min/max
        y_scaleMax(1+cStep:rows-cStep,1+cStep:cols-cStep,cStep2)=max(y_scaleMax(1+cStep:rows-cStep,1+cStep:cols-cStep),tempAvMax);
        y_scaleMin(1+cStep:rows-cStep,1+cStep:cols-cStep,cStep2)=min(y_scaleMin(1+cStep:rows-cStep,1+cStep:cols-cStep),tempAvMin);

        %---- find the derivatives of the max/min envelope at the current level (including all levels so
        %far)

        % find the absolute min/max at all scales

        yMin                                = min(y_scaleMin(:,:,1:1:cStep2),[],3);
        yMax                                = max(y_scaleMax(:,:,1:1:cStep2),[],3);
        %

        sizeFR_S                            = cStep;
        sizeFC_S                            = sizeFR_S;  
        filtG_small                         = gaussF(sizeFR_S,sizeFC_S,1);

        %-----------------------------This is the most time consuming step ------------------------------
        %------------ since this is only a criterion to stop the process before it reaches the ----------
        %------------ last step, subsample to reduce complexity -----------------------------------------
        %---- Low pass filter to adjust the effects of high/low points
        
        yMin0                               = padData(yMin(1:2:end,1:2:end),ceil(sizeFR_S/2));
        yMin1                               = conv2(yMin0,filtG_small);
        yMin2                               = yMin1(sizeFR_S+1:end-sizeFR_S,sizeFC_S+1:end-sizeFC_S);
        yMax0                               = padData(yMax(1:2:end,1:2:end),ceil(sizeFR_S/2));
        yMax1                               = conv2(yMax0,filtG_small);
        yMax2                               = yMax1(sizeFR_S+1:end-sizeFR_S,sizeFC_S+1:end-sizeFC_S);
        
        


        %---- find the derivatives of the max/min envelope
        y_rderiv_Max                        = diff(yMax2,1,1);
        y_cderiv_Max                        = diff(yMax2,1,2);
        y_rderiv_Min                        = diff(yMin2,1,1);
        y_cderiv_Min                        = diff(yMin2,1,2);
        %
        %magnitude of the gradient
        y_magGrad_Max                       = sqrt(y_rderiv_Max(:,2:end).^2+y_cderiv_Max(2:end,:).^2);
        y_magGrad_Min                       = sqrt(y_rderiv_Min(:,2:end).^2+y_cderiv_Min(2:end,:).^2);
        %

        tot_grad_max1(cStep)                = sum(sum(y_magGrad_Max)); %#ok<AGROW>
        tot_grad_min1(cStep)                = sum(sum(y_magGrad_Min)); %#ok<AGROW>

        %
        if (cStep>1)

            diffGradMax                     =  abs((tot_grad_max1(end)-tot_grad_max1(end-2))/tot_grad_max1(end-2));
            diffGradMin                     =  abs((tot_grad_min1(end)-tot_grad_min1(end-2))/tot_grad_min1(end-2));

            %disp([cStep diffGradMin diffGradMax ] )
            if (diffGradMax<stopCriterion)||(diffGradMin<stopCriterion)
                break;
            end
        end
        %
    end

    %
    tot_grad_max                            = tot_grad_max1(end);
    tot_grad_min                            = tot_grad_min1(end);
    %%
    %compare the gradients to decide which surface to keep (smallest) but first re calculate the smooth


    yMin0                               = padData(yMin(1:1:end,1:1:end),ceil(sizeFR_S/2));
    yMin1                               = conv2(yMin0,filtG_small);
    yMin2                               = yMin1(sizeFR_S+1:end-sizeFR_S,sizeFC_S+1:end-sizeFC_S);
    yMax0                               = padData(yMax(1:1:end,1:1:end),ceil(sizeFR_S/2));
    yMax1                               = conv2(yMax0,filtG_small);
    yMax2                               = yMax1(sizeFR_S+1:end-sizeFR_S,sizeFC_S+1:end-sizeFC_S);


    if abs(((tot_grad_max-tot_grad_min)/tot_grad_max))<0.05
        yProm                               = 0.5*yMin2+0.5*yMax2;
        % dataOut(:,:,currLev)=dataIn(:,:,currLev)-(yProm)+mean(yProm(:));
        errSurface                          =  (yProm)-mean(yProm(:));
        dataOut                             = dataIn-errSurface;
    else
        if tot_grad_max>tot_grad_min
            % dataOut(:,:,currLev)=dataIn(:,:,currLev)-(yMin2)+mean(yMin2(:));
            %dataOut=dataIn-(yMin2)+mean(yMin2(:));
            errSurface                      =  (yMin2)-mean(yMin2(:));
            dataOut                         = dataIn-errSurface;
        else
            % dataOut(:,:,currLev)=dataIn(:,:,currLev)-(yMax2)+mean(yMax2(:));
            errSurface                      =  (yMax2)-mean(yMax2(:));
            dataOut                         = dataIn-errSurface;
            %%dataOut=dataIn-(yMax2)+mean(yMax2(:));
        end
    end
    if nargout==4
        errSurfaceMax                       = (yMax2);
        errSurfaceMin                       = (yMin2);
    end
    if reconvertTo255 ==1
        dataOut(dataOut>255)                    = 255;
        dataOut(dataOut<0)                      = 0;
    end
    %dataOut                     = (dataOut-min(dataOut(:)));
    %dataOut                     = 255*(dataOut/max(dataOut(:)));
end


