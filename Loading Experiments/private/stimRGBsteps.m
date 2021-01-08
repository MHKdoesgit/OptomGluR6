

function [ newStimpara ] = stimRGBsteps(desc, stimpara)
%
%%% stimRGBsteps %%%
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
    
    case {'rgbsteps50_contrast1redblack'
            'rgbsteps30contrast50redblack'}
        stimpara.usered = true;
        stimpara.usegreen =false;
        stimpara.useblue = false;
        
    case {'rgbsteps30_contrast50green'
            'rgbsteps30contrast50greenblack'
            'rgbsteps30contrast1greenblack'
            'rgbsteps_greenblack_1greencontrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = false;
        
    case 'rgbsteps_greenblack_0.7greencontrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = false;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0.7;
        
    case {'rgbsteps50_contrast1blueblack'
            'rgbsteps30contrast50blueblack'
            'rgbsteps30contrast1blueblack'
            'rgbsteps_uvblack_1uvcontrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =false;
        stimpara.useblue = true;
        
    case 'rgbsteps_uvblack_0.7uvcontrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =false;
        stimpara.useblue = true;
        stimpara.greencontrast = 0;
        stimpara.bluecontrast = 0.7;
        
    case {'rgbsteps30contrast1cyanblack'
            'rgbsteps_cyanblack_1contrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        
    case 'rgbsteps_cyanblack_0.7contrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0.7;
        
        
    case {'rgbsteps30_contrast50purple'
            'rgbsteps30contrast50purpleblack'}
        stimpara.useblue = true;
        stimpara.usegreen =false;
        stimpara.usered = true;
        
    case {'rgbstepsopponent50_contrast1redblue'
            'rgbsteps30contrast50redblueoppnent' }
        stimpara.opponency = true;
        stimpara.usered = true;
        stimpara.usegreen =false;
        stimpara.useblue = true;
        
    case {'rgbsteps30contrast1bluegreenoppnent'
            'rgbsteps_opponent_1contrast_30frames'}
        stimpara.opponency = true;
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        
    case 'rgbsteps_opponent_0.7contrast_30frames'
        stimpara.opponency = true;
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0.7;
        
    case {'rgbsteps30contrast1bluemeanisoresponse'
            'rgbsteps_uvmean_1uvcontrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0;
        stimpara.bluecontrast = 1;
        stimpara.isoResponse = true;
        
    case 'rgbsteps_uvmean_0.7uvcontrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0;
        stimpara.bluecontrast = 0.7;
        stimpara.isoResponse = true;
        
    case {'rgbsteps30contrast1greenmeanisoresponse'
            'rgbsteps_greenmean_1greencontrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 1;
        stimpara.bluecontrast = 0;
        stimpara.isoResponse = true;
        
    case 'rgbsteps_greenmean_0.7greencontrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0;
        stimpara.isoResponse = true;
        
    case {'rgbsteps30contrast1allrandomized','rgbsteps_allrandomized_1contrast_30frames'}
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 1;
        stimpara.bluecontrast = 1;
        stimpara.seed = -1000;
        
    case 'rgbsteps_allrandomized_0.7contrast_30frames'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0.7;
        stimpara.seed = -1000;
        
    case 'rgbsteps_cyanblack_35contrast_30frames90pfr'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.35;
        stimpara.bluecontrast = 0.35;
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        
    case 'rgbsteps_cyanblack_35contrast_60frames90pfr'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.greencontrast = 0.35;
        stimpara.bluecontrast = 0.35;
        stimpara.nframes = 60;
        stimpara.preframes = 90; 
        
    case 'rgbsteps_greenblack_80contrast_60frames90pfr'
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = false;
        stimpara.greencontrast = 0.8;
        stimpara.bluecontrast = 0;
        stimpara.nframes = 60;
        stimpara.preframes = 90; 
        stimpara.coneisolating = false;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        
        
    case 'rgbsteps_uvblack_80contrast_60frames90pfr'
        stimpara.usered = false;
        stimpara.usegreen =false;
        stimpara.useblue = ture;
        stimpara.greencontrast = 0;
        stimpara.bluecontrast = 0.8;
        stimpara.nframes = 60;
        stimpara.preframes = 90; 
        stimpara.coneisolating = false;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
                
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end


