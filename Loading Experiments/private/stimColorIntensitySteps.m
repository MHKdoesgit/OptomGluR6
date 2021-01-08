

function [ newStimpara ] = stimColorIntensitySteps(desc, stimpara)
%
%%% stimColorIntensitySteps %%%
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
    
    % Color Intensity Steps
    case 'colorintensitysteps120bluebkg30redsteps8skipval'        
        stimpara.redmax = false;
        stimpara.bluemax = true;
        stimpara.framesonscreen = 30;
        stimpara.skipvalue = 8;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorintensitysteps120redbkg30bluesteps8skipval'        
        stimpara.redmax = true;
        stimpara.bluemax = false;
        stimpara.framesonscreen = 30;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 8;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorintensitysteps60bluebkg30redsteps8skipval'        
        stimpara.redmax = false;
        stimpara.bluemax = true;
        stimpara.framesonscreen = 30;
        stimpara.backgroundonscreen = 60;
        stimpara.skipvalue = 8;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorintensitysteps60redbkg30bluesteps8skipval'        
        stimpara.redmax = true;
        stimpara.bluemax = false;
        stimpara.framesonscreen = 30;
        stimpara.backgroundonscreen = 60;
        stimpara.skipvalue = 8;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
    case {'colorintensitysteps120uvbkg60greensteps16skipval'
            'colorintensitysteps_120uvbkg_60greensteps_16skipval'}        
        stimpara.greenmax = false;
        stimpara.bluemax = true;
        stimpara.framesonscreen = 60;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 16;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case {'colorintensitysteps120greenbkg60uvsteps16skipval'
            'colorintensitysteps120greebbkg60uvsteps16skipval'
            'colorintensitysteps_120greenbkg_60uvsteps_16skipval'},        
        stimpara.greenmax = true;
        stimpara.bluemax = false;
        stimpara.framesonscreen = 60;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 16;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorintensitysteps_120uvbkg_60greensteps_10skipval'        
        stimpara.greenmax = false;
        stimpara.bluemax = true;
        stimpara.framesonscreen = 60;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 10;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.maxcontrsat = 0.7;
        
    case 'colorintensitysteps_120greenbkg_60uvsteps_10skipval'        
        stimpara.greenmax = true;
        stimpara.bluemax = false;
        stimpara.framesonscreen = 60;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 10;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.maxcontrsat = 0.7;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end


