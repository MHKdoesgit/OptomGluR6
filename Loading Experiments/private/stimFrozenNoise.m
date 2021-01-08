

function [ newStimpara ] = stimFrozenNoise (desc, stimpara)
%
%%% stimFrozenNoise %%%
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
    
    case 'frozennoise1600x1200gauss_rn1500_fn3002blinks'
        
        colormode = questdlg('Choose the stimulus color','Black & white or Color','Black & white','Color','Color');
        switch colormode
            case 'Color'
                stimpara.RunningFrames = 1500;
                stimpara.FrozenFrames = 300;
                stimpara.color = true;
                stimpara.independentcolors = true;
                [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
                [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.3,0.3);
                stimpara.coneisolating = true;
                stimpara.seedrunningred = -1000;
                stimpara.seedrunninggreen = -1000/2;
                stimpara.seedrunningblue = -1000/4;
                stimpara.seedfrozenred = -10000;
                stimpara.seedfrozengreen = -10000/2;
                stimpara.seedfrozenblue = -10000/4;
                
            case 'Black & white'
                stimpara.nblinks = 2;
                stimpara = rmfield(stimpara,{'color'});
                
        end
        
    case 'frozennoise_1600x1200_gauss_4blinks_1200run900frozen'
        stimpara.RunningFrames = 1200;
        stimpara.FrozenFrames = 900;
        stimpara.blackwhite = false;
        stimpara.nblinks = 4;
        stimpara = rmfield(stimpara,{'color'});
        
        % single color used for OptomGluR6 project
    case {'frozennoise1600x1200gauss_rn1500_fn3001blinks','opto_frozennoise1600x1200gauss_rn1500_fn3001blinks'}
        stimpara.nblinks = 1;
        stimpara = rmfield(stimpara,{'color'});
        
    case {'opto_frozennoise1600x1200gauss_rn1500_fn3002blinks'}
        stimpara.nblinks = 2;
        stimpara = rmfield(stimpara,{'color'});
        
    case 'checkerflicker_2x2bw4blinks_4200run_300frozen'
        stimpara.RunningFrames = 4200;
        stimpara.FrozenFrames = 300;
        stimpara.blackwhite = false;
        stimpara.nblinks = 4;
        stimpara = rmfield(stimpara,{'color'});
        
    case 'frozennoise1600x1200gauss_rn1500_fn300_4blinks'
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 4;
        stimpara.color = true;
        stimpara.independentcolors = true;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.3,0.3);
        stimpara.coneisolating = true;
        stimpara.seedrunningred = -1000;
        stimpara.seedrunninggreen = -1000/2;
        stimpara.seedrunningblue = -1000/4;
        stimpara.seedfrozenred = -10000;
        stimpara.seedfrozengreen = -10000/2;
        stimpara.seedfrozenblue = -10000/4;
        
    case 'frozennoise1600x1200gauss_rn1500_fn300_2blinks'
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 2;
        stimpara.color = true;
        stimpara.independentcolors = true;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.3,0.3);
        stimpara.coneisolating = true;
        stimpara.seedrunningred = -1000;
        stimpara.seedrunninggreen = -1000/2;
        stimpara.seedrunningblue = -1000/4;
        stimpara.seedfrozenred = -10000;
        stimpara.seedfrozengreen = -10000/2;
        stimpara.seedfrozenblue = -10000/4;
        
    case {'frozennoise1600x1200gauss_rn1500_fn300_1blink','frozennoise1600x1200gauss_rn1500_fn300_1blinks',...
            'frozennoise_1600x1200_gauss_rn1500_fn300_color_1blink'}
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 1;
        stimpara.color = true;
        stimpara.independentcolors = true;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0,0.3,0.3);
        stimpara.coneisolating = true;
        stimpara.seedrunningred = -1000;
        stimpara.seedrunninggreen = -1000/2;
        stimpara.seedrunningblue = -1000/4;
        stimpara.seedfrozenred = -10000;
        stimpara.seedfrozengreen = -10000/2;
        stimpara.seedfrozenblue = -10000/4;
        
    case 'frozennoise1600x1200gauss_rgb_rn1500_fn300_1blink'
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 1;
        stimpara.color = true;
        stimpara.independentcolors = true;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0.5,0.5,0.5);
        [stimpara.redcontrast,stimpara.greencontrast ,stimpara.bluecontrast] = deal(0.3,0.3,0.3);
        stimpara.coneisolating = false;
        stimpara.seedrunningred = -1000;
        stimpara.seedrunninggreen = -1000/2;
        stimpara.seedrunningblue = -1000/4;
        stimpara.seedfrozenred = -10000;
        stimpara.seedfrozengreen = -10000/2;
        stimpara.seedfrozenblue = -10000/4;
        
    case 'fffgausswithrepeats_2blinks'
        stimpara.seedrunningnoise = -10000;
        stimpara.seedfrozennoise = -20000;
        stimpara.RunningFrames = 900;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 2;
        stimpara.pulserate = 2;
        stimpara.blackwhite = false;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.3;
        stimpara = rmfield(stimpara,{'color'});
        
        %========================checker flicker with frozen and running noise ====================%
        
    case 'frozennoise8x8bw1blink1500run300freeze'
        
        colormode = questdlg('Choose the stimulus color','Black & white or Color','Black & white','Color','Color');
        switch colormode
            case 'Color'
                stimpara.stixelwidth = 8;
                stimpara.stixelheight = 8;
                stimpara.RunningFrames = 1500;
                stimpara.FrozenFrames = 300;
                stimpara.nblinks = 1;
                stimpara.color = true;
                stimpara.independentcolors = true;
                [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
                [stimpara.redContrast,stimpara.greenContrast ,stimpara.blueContrast] = deal(0,0.7,0.7);
                [stimpara.usered,stimpara.usegreen ,stimpara.useblue] = deal(0,1,1);
                stimpara.coneisolating = true;
                stimpara.seedrunningred = -1000;
                stimpara.seedrunninggreen = -1000/2;
                stimpara.seedrunningblue = -1000/4;
                stimpara.seedfrozenred = -10000;
                stimpara.seedfrozengreen = -10000/2;
                stimpara.seedfrozenblue = -10000/4;
            case 'Black & white'
                stimpara.color = false;
                stimpara.stixelwidth = 8;
                stimpara.stixelheight = 8;
                stimpara.RunningFrames = 1500;
                stimpara.FrozenFrames = 300;
                stimpara.nblinks = 1;
                stimpara.blackwhite = true;
                stimpara.meanintensity = 0.5;
                stimpara.contrast = 1;
                stimpara = rmfield(stimpara,{'color'});
        end
        
    case {'frozennoise8x8binary_rn1500_fn300_1blinks'}
        stimpara.color = false;
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 1;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara = rmfield(stimpara,{'color'});
        
    case 'frozennoise_8x8_bw_rn1500_fn300_color_1blink'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 1;
        stimpara.color = true;
        stimpara.independentcolors = true;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        [stimpara.redContrast,stimpara.greenContrast ,stimpara.blueContrast] = deal(0,0.7,0.7);
        [stimpara.usered,stimpara.usegreen ,stimpara.useblue] = deal(0,1,1);
        stimpara.coneisolating = true;
        stimpara.seedrunningred = -1000;
        stimpara.seedrunninggreen = -1000/2;
        stimpara.seedrunningblue = -1000/4;
        stimpara.seedfrozenred = -10000;
        stimpara.seedfrozengreen = -10000/2;
        stimpara.seedfrozenblue = -10000/4;        
        
    case {'opto_frozennoise8x8bw1blink1500run300freeze'}
        stimpara.color = false;
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 1;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara = rmfield(stimpara,{'color'});
        
    case {'frozennoise8x8bw2blinks1500run300freeze','opto_frozennoise8x8bw2blinks1500run300freeze'}
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 2;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        
    case 'frozennoise8x8bw1blink3000run600freeze'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.RunningFrames = 3000;
        stimpara.FrozenFrames = 600;
        stimpara.nblinks = 1;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        
    case {'frozennoise5x5bw1blink1500run300freeze','frozennoisex5bw1blink1500run300freeze',...
            'frozennoise5x5bwrgb1blink1500run300freeze','frozennoise5x5bw1blink_1500run_300freeze',...
            'frozennoise5x5bw1blink1500run300freeze_after','frozennoise5x5bw1blink1500runs300freeze',...
            'frozennoise5x5bw1blinkrgb1500run300freeze'}
        
        colormode = questdlg('Choose the stimulus color','Black & white or Color','Black & white','Color','Color');
        switch colormode
            case 'Color'
                stimpara.stixelwidth = 5;
                stimpara.stixelheight = 5;
                stimpara.RunningFrames = 1500;
                stimpara.FrozenFrames = 300;
                stimpara.nblinks = 1;
                stimpara.color = true;
                stimpara.independentcolors = true;
                [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
                [stimpara.redContrast,stimpara.greenContrast ,stimpara.blueContrast] = deal(0,0.7,0.7);
                [stimpara.usered,stimpara.usegreen ,stimpara.useblue] = deal(0,1,1);
                stimpara.coneisolating = true;
                stimpara.seedrunningred = -1000;
                stimpara.seedrunninggreen = -1000/2;
                stimpara.seedrunningblue = -1000/4;
                stimpara.seedfrozenred = -10000;
                stimpara.seedfrozengreen = -10000/2;
                stimpara.seedfrozenblue = -10000/4;
            case 'Black & white'
                stimpara.color = false;
                stimpara.stixelwidth = 5;
                stimpara.stixelheight = 5;
                stimpara.RunningFrames = 1500;
                stimpara.FrozenFrames = 300;
                stimpara.nblinks = 1;
                stimpara.blackwhite = true;
                stimpara.meanintensity = 0.5;
                stimpara.contrast = 1;
                stimpara = rmfield(stimpara,{'color'});
        end
        
    case 'checkerflicker9x9bw_3blinks_withrepeats'
        stimpara.fps = 75;
        stimpara.stixelsize = 2.5;
        stimpara.seedrunningnoise = -10000;
        stimpara.seedfrozennoise = -20000;
        stimpara.color = false;
        stimpara.stixelwidth = 9;
        stimpara.stixelheight = 9;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.nblinks = 3;
        stimpara.pulserate = 3;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara = rmfield(stimpara,{'color'});
        
    case {'frozennoise5x5bw1blink3600run1200freeze'
            'frozennoise5x5bw1blink3600run1200freeze_40br'
            'frozennoise5x5bw1blink3600run1200freeze_0br'}
        
        stimpara.color = false;
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 5;
        stimpara.RunningFrames = 3600;
        stimpara.FrozenFrames = 1200;
        stimpara.nblinks = 1;
        stimpara.blackwhite = true;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara = rmfield(stimpara,{'color'});
        
    case {'fullfieldfrozennoise_gauss_1blink_1800run_600freeze'
            'fullfieldfrozennoise_gauss_1blink_1800run_600freeze_40'
            'fullfieldfrozennoise_gauss_1blink_1800run_600freeze_0'}
        stimpara.color = false;
        stimpara.stixelwidth = 1600;
        stimpara.stixelheight = 1200;
        stimpara.RunningFrames = 1800;
        stimpara.FrozenFrames = 600;
        stimpara.nblinks = 1;
        stimpara.blackwhite = false;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.3;
        stimpara = rmfield(stimpara,{'color'});  
end

newStimpara = stimpara;

end

