

function [ newStimpara ] = stimCheckerflicker(desc,stimpara)
%
%%% stimCheckerflicker %%%
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

if (stimpara.blackwhite)
    stimpara.contrast = 1;
else
    stimpara.contrast = 0.3;
end

switch lower(desc)
    
    case {'checkerflicker10x10bw2blinks'
            'checkerflicker10x10bw2blinks2'
            'checkerflicker10x0bw_2blinks'
            'checkerflicker10x10bw_2blinks'
            'checkerflickerbw10x10_2bl'
            'checkerflicker_10x10_bw_2blinks'
            'checkerflicker10x10_bw2blinks'
            'checkerflicker_10x10bw_2blinks'
            'checkerflickergabazine_10x10_bw_2blinks'
            'checkerflickertpmpa_10x10_bw_2blinks'
            'checkerflickerlap4_10x10_bw_2blinks'
            'checkerflickerstrychnine_10x10_bw_2blinks'
            'checkerflickerwashout_10x10_bw_2blinks'}
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 10;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case {'checkerflicker8x8_bw2blinks'
            'checkerflicker8x8bw2blincks'
            'checkerflicker8x8bw2blinksblackwhite'
            'checkerflicker_8x8_bw_2blinks'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case {'checkerflicker_8x8_bw_1blinks','checkerflicker_8x8_bw_1contrast_1blink',...
            'opto_checkerflicker_8x8_bw_1contrast_1blink'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.nblinks = 1;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case {'checkerflicker5x5bw2blinks'
            'checkerflicker5x5_2blinks'}
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 5;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case {'checkerflicker5x51blink'}
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 5;
        stimpara.nblinks = 1;
        stimpara.seed = -10000;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker1x1bw2blinks'
        stimpara.stixelwidth = 1;
        stimpara.stixelheight = 1;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker_8x8_bw_4blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.nblinks = 4;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker_3x3_bw_6blinks'
        stimpara.stixelwidth = 3;
        stimpara.stixelheight = 3;
        stimpara.nblinks = 6;
        stimpara.seed = -1000;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker_2x2bw_4blinks'
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 4;
        stimpara.seed = -10000;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker_3x3bw_4blinks'
        stimpara.stixelwidth = 3;
        stimpara.stixelheight = 3;
        stimpara.nblinks = 4;
        stimpara.seed = -10000;
        stimpara = rmfield(stimpara,{'color','coneisolating','secondseed','thirdseed'});
        
    case 'checkerflicker10x10bw2blinkscolor'
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 10;
        stimpara.color= true;
        stimpara.usered = true;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 1;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        
    case {'checkerflicker8x8bw2blinks color'
            'checkerflicker8x8bw2blinkscolor'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        
    case {'checkerflicker6x6bw2blinkscolor'
            'checkerflicker6x6bw2blinkscolor2seed'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        
    case 'checkerflicker6x6gauss2blinkscolor'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.blackwhite = false;
        
        
    case 'checkerflicker3x3bw2blinkscolor'
        stimpara.stixelwidth = 3;
        stimpara.stixelheight = 3;
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
        
    case 'checkerflicker2x2bw2blinkscolor'
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
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
        
    case 'checkerflicker1x1bw2blinkscolor'
        stimpara.stixelwidth = 1;
        stimpara.stixelheight = 1;
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
        
    case 'checkerflicker8x8bw2blinksredonly'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.color= true;
        stimpara.usered = true;
        stimpara.usegreen = false;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 1;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0;
        
    case 'checkerflicker8x8bw2blinksgreenonly'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        
    case 'checkerflicker6x6bw2blinksgreenonly'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        
    case {'checkerflicker8x8bw2blinksuvonly'
            'checkerflicker8x8bw2blinksblueonly'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        
    case {'checkerflicker6x6bw2blinksblueonly'
            'checkerflicker6x6bw2blinksuvonly'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 2;
        
        
    case 'checkerflicker8x8bw2blinkscyanonly'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        
        
    case 'checkerflicker6x6bw2blinkscyanonly'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        
        %%%% new format
    case 'checkerflicker_2x2_bw_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_2x2_bw_color_0.7contrast_4blinks'
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 4;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_4x4_bw_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_4x4_bw_cyanonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_4x4_bw_greenonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_4x4_bw_uvonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_4x4_gauss_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
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
        stimpara.coneisolating = true;
        stimpara.blackwhite = false;
        
    case 'checkerflicker_6x6_bw_color_0.7contrast_1blink'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.coneisolating = true;
        stimpara.nblinks = 1;
        
    case 'checkerflicker_6x6_bw_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_6x6_bw_color_0.7contrast_6blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.coneisolating = true;
        stimpara.nblinks = 6;
        
    case 'checkerflicker_6x6_bw_cyanonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_6x6_bw_greenonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_6x6_bw_uvonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_6x6_gauss_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
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
        stimpara.coneisolating = true;
        stimpara.blackwhite = false;
        
    case 'checkerflicker_8x8_bw_color_0.7contrast_2blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_8x8_bw_color_0.7contrast_1blink'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        stimpara.coneisolating = true;
        stimpara.nblinks = 1;
        
    case 'checkerflicker_4x4_bw_color_0.7contrast_2blink'
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
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
        stimpara.coneisolating = true;
        stimpara.nblinks = 2;
        
    case 'checkerflicker_8x8_bw_cyanonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        stimpara.coneisolating = true;
        
    case 'checkerflicker_8x8_bw_greenonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = false;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0.7;
        stimpara.blueContrast = 0;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_8x8_bw_uvonly_0.7contrast_2blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.color= true;
        stimpara.usered = false;
        stimpara.usegreen = false;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 0;
        stimpara.greenContrast = 0;
        stimpara.blueContrast = 0.7;
        stimpara.coneisolating = true;
        
    case 'checkerflicker_8x8_gauss_color_0.3contrast_2blinks'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
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
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.coneisolating = true;
        stimpara.blackwhite = false;
        
    case 'checkerflicker_5x5_bw_3colors_1contrast_1blink'
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 5;
        stimpara.color= true;
        stimpara.usered = true;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0.5;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redContrast = 1;
        stimpara.greenContrast = 1;
        stimpara.blueContrast = 1;
        stimpara.seed = -1000;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.coneisolating = false;
        stimpara.blackwhite = true;
        stimpara.nblinks = 1;
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;
end

