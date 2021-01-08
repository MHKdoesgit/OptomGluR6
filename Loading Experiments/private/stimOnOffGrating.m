

function [ newStimpara ] = stimOnOffGrating (desc, stimpara)
%
%%% stimOnOffGrating %%%
%
%
% This function is a private function of loadStimulusParameters function.
%
%================================Inputs====================================
%
%   desc : original name of the mcd file and experiment data.
%   stimpara : stimulus parameter structure.
%
%================================Output====================================
%
%   newStimpara : updated stimpara structure.
%

switch lower(desc)
    
    case {'onoffgrating60_for60hz'
            'onoffgrating60'
            'onoffgrating60d'
            'onoff30grating60d'}
        stimpara.meanintensity  = 0.5;
        stimpara.contrast  = 1;
        
    case {'onoffgratingp60s30ph30ang0bluegreen'}
        stimpara.period = 60;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.7,-0.7);
        stimpara.coneisolating = true;
        
    case {'onoffgratingp60s30ph30ang0cyanblack'}
        stimpara.period = 60;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.7,0.7);
        stimpara.coneisolating = true;
        
    case {'onoffgratingp60s30ph30ang0greenblack'}
        stimpara.period = 60;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.7,0);
        stimpara.coneisolating = true;
        
    case {'onoffgratingp60s30ph30ang0blueblack'}
        stimpara.period = 60;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.0,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0,0.7);
        stimpara.coneisolating = true;
        
    case {'onoffgratingp32s30ph30ang0cyanblack'}
        stimpara.period = 32;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.7,0.7);
        stimpara.coneisolating = true;
end

newStimpara = stimpara;

end


