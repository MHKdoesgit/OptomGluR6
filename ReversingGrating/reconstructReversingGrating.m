function [ stimRecon ] = reconstructReversingGrating( stimPara, nx, ny)
%RECONSTRUCTREVERSINGGRATING Summary of this function goes here
%   Detailed explanation goes here
%have to fix what happens to the end of the screen in the reversing phase
%have to find out what is happening with rounding

phases      =       stimPara.Nphases;
widths      =       stimPara.stripewidths;

stimRecon   =       stimPara.meanintensity*ones(sum(phases), ny, nx);
fact        =       stimPara.contrast*stimPara.meanintensity;
stimRecon   =       stimRecon-fact;

idx         =       1;

for iWidth  =   1:numel(widths)
    
    whiteStarts     =   1:2*widths(iWidth):nx;
    myMat           =   zeros( ny, nx);
    
    for ii  =   1:numel(whiteStarts)
        myMat(:,whiteStarts(ii):whiteStarts(ii)+widths(iWidth)-1) = stimPara.meanintensity + fact;
    end
    nPhases =       phases(iWidth);
    
    for iPhase      =   0:nPhases-1
        newMat      =   [myMat(:,round(iPhase*widths(iWidth)/nPhases+1:end))...
            (stimPara.meanintensity-fact)*ones(ny,round(iPhase*widths(iWidth)/nPhases))];
        stimRecon(idx,:,:)      =       newMat(:,1:nx);
        idx         =   idx+1;
    end
end

end

