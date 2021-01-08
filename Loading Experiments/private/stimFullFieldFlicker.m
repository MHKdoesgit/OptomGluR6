

function [ newStimpara ] = stimFullFieldFlicker(desc,stimpara)
%
%%% stimFullFieldFlicker %%%
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

% correct the contrast
if (stimpara.blackwhite)
    stimpara.contrast = 1;
else
    stimpara.contrast = 0.3;
end

switch lower(desc)
    
    case {'fff2blinks'
            'fff2bl'
            'fff_gauss_2blinks'
            'fff_gauss_2blinks_2nd'
            'fff_gauss_2blinks_3rd'
            'ff_checkerflicker1600x1200_ga_2blinks'
            'checkerflicker1600x1200gauss2blinks'
            'checkerflicker1600x1200_ga_2blinks'}
        stimpara.nblinks = 2;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating'});
        
    case {'fff1blink', 'fff1blinks','fff_gauss1blink','fff_gauss1blink_bright79'...
            'fff1blink_bright91','fffgauss1blink','fff1blinkafter','fff1blink_duringadapt'}
        stimpara.nblinks = 1;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating'});
        if strcmpi(desc,'fff_gauss1blink_bright79')
            stimpara.OLEDbrightness = 79;
        elseif strcmpi(desc,'fff1blink_bright91')
            stimpara.OLEDbrightness = 91;
        end
        
    case 'fullfiledflicker_1600x1200_gauss_4blinks'
        stimpara.nblinks = 4;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating'});
        
    case {'fullfieldflicker_1600x1200_gauss_1blink'
            'opto_fullfieldflicker_1600x1200_gauss_1blink'}
        stimpara.nblinks = 1;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating'});
        
    case {'checkerflicker1600x1200gauss2blinkscolor'
            'fullfieldflicker_1600x1200_gauss_color_2blinks'}
        stimpara.nblinks = 2;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.3;
        stimpara.blueContrast = 0.3;
        stimpara.seed = -1000;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.coneisolating = true;
        
    case{'checkerflicker1600x1200gauss1blinks color'
            'checkerflicker1600x1200gauss1blinkcolor'
            'checkerflicker1600x1200gauss1blinkscolor'
            'fullfieldflicker_1600x1200_gauss_color_1blinks'}
        stimpara.nblinks = 1;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.3;
        stimpara.blueContrast = 0.3;
        stimpara.seed = -1000;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.coneisolating = true;
        
    case {'checkerflicker1600x1200gauss2blinkscyanonly'
            'fullfieldflicker_1600x1200_gauss_cyanonly_2blinks'}
        stimpara.nblinks = 2;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.3;
        stimpara.blueContrast = 0.3;
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        
    case 'checkerflicker1600x1200gauss2blinksredonly'
        stimpara.nblinks = 2;
        stimpara.color= true;
        stimpara.usered = true;
        stimpara.usegreen = false;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0.3;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0;
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        
    case 'checkerflicker1600x1200gauss1blinksredonly'
        stimpara.nblinks = 1;
        stimpara.color= true;
        stimpara.usered = true;
        stimpara.usegreen = false;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0.3;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0;
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        
    case {'checkerflicker1600x1200gauss2blinksgreenonly'
            'fullfieldflicker_1600x1200_gauss_greenonly_2blinks'}
        stimpara.nblinks = 2;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.3;
        stimpara.blueContrast = 0;
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        
    case {'checkerflicker1600x1200gauss2blinksuvonly'
            'checkerflicker1600x1200gauss2blinksblueonly'
            'fullfieldflicker_1600x1200_gauss_uvonly_2blinks'}
        stimpara.nblinks = 2;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.3;
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        
    case 'checkerflicker1600x1200gauss1blinksblueonly'
        
        stimpara.nblinks = 1;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end


