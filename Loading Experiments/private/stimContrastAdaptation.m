

function [ newStimpara ] = stimContrastAdaptation (desc, stimpara)
%
%%% stimContrastAdaptation %%%
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
    
    case {'ff_ca_ga_32_08_1blinks_sf2400'
            'ff_ca_ga_32_08_1blink_2400frames'
            'ff_ca_ga_32_08_1blinks2400frame'
            'ff_ca_ga_32_08_1blinks2400frames'
            'ff_ca_ga_32_08_1blinks_2400sf'},
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.lowcontrast = 0.08;
        
    case {'ff_ca_ga_32_08_2blinks_sf1200'
            'ff_ca_ga_32_08_2blinks'
            'ff_ca_ga_32_08_2blink_2400frames'},
        stimpara.lowcontrast = 0.08;
        
    case {'ff_ca_ga_32_04_1blinks_sf2400'
            'ff_ca_ga_32_04_1blink_2400frames'
            'ff_ca_ga_32_04_1blinks_2400frames'
            'ff_ca_ga_32_04_1blinks2400frames'},
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.lowcontrast = 0.04;
        
    case {'ff_ca_ga_32_04_2blinks_sf1200'
            'ff_ca_ga_32_04_2blink_2400frames'
            'ff_ca_ga_32_04_2blinks'},
        stimpara.lowcontrast = 0.04;
        
    case {'ff_ca_ga_32_16_1blinks_sf2400'
            'ff_ca_ga_32_16_1blink_2400frames'
            'ff_ca_ga_32_16_1blinks_2400frames'
            'ff_ca_ga_32_16_1blinks2400frames'},
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.lowcontrast = 0.16;
        
    case {'ff_ca_ga_32_16_2blinks_sf1200'
            'ff_ca_ga_32_16_2blinks'
            'ff_ca_ga_32_16_2blink_2400frames'},
        stimpara.lowcontrast = 0.16;
        
    case 'ff_ca_32_12_2blinks_sf1200_2blinks'
        stimpara.lowcontrast = 0.12;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc]);
end;

newStimpara = stimpara;

end


