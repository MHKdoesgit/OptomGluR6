

function [ newStimpara ] = stimRotatingStrips(desc, stimpara)
%
%%% stimRotatingStrips %%%
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
    
    case 'rotatingstrips6x6bw2blinkscolor2seed'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 1;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case {'rotatingstrips6x6bw2blinkscyanonly'
            'rotatingstrips6x6bw2blinks2seedcyanonly'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 1;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case {'rotatingstrips6x6bw2blinksgreenonly'
            'rotatingstrips6x6bw2blinks2seedgreenonly'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 0;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case {'rotatingstrips6x6bw2blinksuvonly'
            'rotatingstrips6x6bw2blinks2seeduvonly'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 1;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case {'rotating_stripes_2x600_2blinks'
            'rotatingstripes_2x600_bw_120preframes_18000dur'}
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;        
        
    case {'rotating_stripes_4x600_2blinks'
            'rotatingstripes_4x600_bw_120preframes_18000dur'}
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        
    case {'rotatingstripes5bw1blink5min'}
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        stimpara.nblinks = 1;
        
    case {'rotatingstripes_4x600_1blinks','rotating_stripes_4x600_1blinks','rotatingstripes4x600_1blink'}
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        stimpara.seed = -1000;
        stimpara.nblinks = 1;
        
    case {'rotatingstripes_2x600_1blinks','rotatingstripes2x600_1blinks','rotating_stripies_2x600_1blinks'...
            ,'rotatingstripes_2x600_1blink','rotating_stripes_2x600_1blinks'}
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        stimpara.seed = -1000;
        stimpara.nblinks = 1;
        
    case 'rotatingstripes_4x600_4blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        stimpara.seed = -1000;
        
    case {'rotating_stripes_8x600_2blinks'
            'rotatingstripes_8x600_bw_120preframes_18000dur'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.stimduration = 18000;
        stimpara.preframes = 120;
        
    case 'rotatingstrips6x6gauss2blinkscolor2seed'
        stimpara.blackwhite = false;
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 1;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_color_120preframes_18000dur'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_color_120pfr_18000dur4blinks'
        stimpara.nblinks = 4;
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_color_120preframes_36000dur'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 36000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_color_120pfr_36000dur4blinks'
        stimpara.nblinks = 4;
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 36000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_cyanonly_120preframes_18000dur'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_4x600_bw_greenonly_120preframes_18000dur'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case    'rotatingstripes_4x600_bw_uvonly_120preframes_18000dur'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_6x600_bw_color_120preframes_18000dur'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_6x600_bw_cyanonly_120preframes_18000dur'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_6x600_bw_greenonly_120preframes_18000dur'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case    'rotatingstripes_6x600_bw_uvonly_120preframes_18000dur'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_8x600_bw_color_120preframes_18000dur'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_8x600_bw_cyanonly_120preframes_18000dur'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_8x600_bw_greenonly_120preframes_18000dur'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    case 'rotatingstripes_8x600_bw_uvonly_120preframes_18000dur'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.stimduration = 18000;
        stimpara.coneisolating = true;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end

