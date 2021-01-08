

function [ newStimpara ] = stimDirectionGrating (desc, stimpara)
%
%%% stimDirectionGrating %%%
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
    case {'directiongratingslow'
            'directiongrating_slow'}
        stimpara.period = 100;
        stimpara.duration = 500;
        stimpara.regeneration = 200;
        stimpara.gratingwidth = 350;
        stimpara.brightness = 1;
        
    case {'directiongratingfast'
            'directiongrating_fast'}
        stimpara.period = 50;
        stimpara.duration = 500;
        stimpara.regeneration = 200;
        stimpara.gratingwidth = 350;
        stimpara.brightness = 1;
        
    case {'direction_grating_bar40_t80_10x5'
            'directiongrating_bar40_t80_10x5'}
        stimpara.period = 80;
        stimpara.regeneration = 200;
        stimpara.gratingwidth = 300;
        stimpara.cycles = 10;
        stimpara.brightness = 1;
        
    case 'direction_grating_fast_2hz_for_60hz'
        stimpara.cycles = 3;
        stimpara.period = 30;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 350;
        
    case 'direction_grating_slow_1hz_for_60hz'
        stimpara.cycles = 3;
        stimpara.period = 60;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 350;
        
    case 'directiongrating_cyc10_bar300_dur400by100_bluegreen'
        stimpara.color = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.coneisolating = true;
        stimpara.colorfirst = [0,0.7,0];
        stimpara.colorsecond = [0,0,0.7];
        
    case 'directiongrating_cyc10_bar300_dur400by100_cyanblack'
        stimpara.color = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.coneisolating = true;
        stimpara.colorfirst = [0,0,0];
        stimpara.colorsecond = [0,0.7,0.7];
        
    case 'directiongrating_cyc10_bar300_dur400by100_bluebalck'
        stimpara.color = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.coneisolating = false;
        stimpara.colorfirst = [0,0,0];
        stimpara.colorsecond = [0,0,0.7];
        
    case 'directiongrating_cyc10_bar300_dur400by100_greenblack'
        stimpara.color = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.coneisolating = false;
        stimpara.colorfirst = [0,0,0];
        stimpara.colorsecond = [0,0.7,0];
        
    case 'directiongrating30_4hz'
        stimpara.cycles = 3;
        stimpara.period = 15;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 112.5;
        
    case 'directiongrating60_4hz'
        stimpara.cycles = 3;
        stimpara.period = 15;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 225;
        
    case  'direction_grating32_4hz'
        stimpara.cycles = 5;
        stimpara.period = 15;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 120;
        
    case 'direction_grating60_2hz'
        stimpara.cycles = 5;
        stimpara.period = 30;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 225;
        
    case 'direction_grating60_4hz'
        stimpara.cycles = 5;
        stimpara.period = 15;
        stimpara.duration = 300;
        stimpara.regeneration = 120;
        stimpara.gratingwidth = 225;
        
    case 'directiongratingsequence_bw_6stimuli'
        stimpara.numstimuli = 6;
        stimpara.cycles = [5, 5, 5, 5, 5, 5];
        stimpara.period = [100, 60, 60, 150, 70, 15];
        stimpara.duration = [200, 240, 240, 600, 280, 180];
        stimpara.regeneration = [300, 120, 120, 300, 120, 120];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1, 1, 0, 0, 1];
        stimpara.gratingwidthwhite = [117.5, 300, 425, 120, 141, 112.5];
        stimpara.gratingwidthblack = [117.5, 900, 425, 120, 141, 112.5];
        stimpara.contrast = [1, 1, 1, 1, 1, 1];
        
    case 'directiongratingsequence_color_4stimuli'
        stimpara.numstimuli = 4;
        stimpara.cycles = [5, 5, 5, 5];
        stimpara.period = [100, 150, 70, 15];
        stimpara.duration = [400, 600, 280, 180];
        stimpara.regeneration = [300, 300, 120, 120];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 0, 0, 1];
        stimpara.gratingwidthwhite = [117.5, 120, 141, 112.5];
        stimpara.gratingwidthblack = [117.5, 120, 141, 112.5];
        stimpara.contrast = [0.7, 0.7, 0.7, 0.7];
        stimpara.meanintensity = [0.0, 0.5, 0.5];
        stimpara.coneIsolating = true;
        
    case 'directiongratingsequence_color_3stimuli'
        stimpara.numstimuli = 3;
        stimpara.cycles = [4, 4, 4];
        stimpara.period = [100, 15, 120];
        stimpara.duration = [400, 180, 480];
        stimpara.regeneration = [300, 120, 120];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1, 1];
        stimpara.gratingwidthwhite = [117.5, 112.5, 480];
        stimpara.gratingwidthblack = [117.5, 112.5, 480];
        stimpara.contrast = [0.7, 0.7, 0.7];
        stimpara.meanintensity = [0.0, 0.5, 0.5];
        stimpara.coneIsolating = true;
        
    case 'directiongratingsequence_bw_3stimuli_1contrast'
        stimpara.numstimuli = 3;
        stimpara.cycles = [4, 4, 4];
        stimpara.period = [100, 15, 120];
        stimpara.duration = [400, 180, 480];
        stimpara.regeneration = [300, 120, 180];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1, 1];
        stimpara.gratingwidthwhite = [117.5, 112.5, 480];
        stimpara.gratingwidthblack = [117.5, 112.5, 480];
        stimpara.contrast = [1, 1, 1];
        stimpara.meanintensity = [0.5, 0.5, 0.5];
        stimpara.coneIsolating = false;
        stimpara.color = false;
        
    case 'directiongratingsequence_color_3stimuli_12angles'
        stimpara.Nangles = 12;
        stimpara.numstimuli = 3;
        stimpara.cycles = [4, 4, 4];
        stimpara.period = [100, 15, 30];
        stimpara.duration = [400, 180, 300];
        stimpara.regeneration = [300, 120, 300];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1, 0];
        stimpara.gratingwidthwhite = [117.5, 112.5, 150];
        stimpara.gratingwidthblack = [117.5, 112.5, 150];
        stimpara.contrast = [1, 1, 1];
        stimpara.meanintensity = [0.5, 0.5, 0.5];
        stimpara.coneIsolating = false;
        stimpara.color = false;
        
    case 'directiongratingsequence_color_3stimuli_8812angles'
        stimpara.Nangles = [8,8,12];
        stimpara.numstimuli = 3;
        stimpara.cycles = [4, 4, 4];
        stimpara.period = [100, 15, 15];
        stimpara.duration = [400, 180, 300];
        stimpara.regeneration = [300, 120, 300];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1, 0];
        stimpara.gratingwidthwhite = [125, 120, 160];
        stimpara.gratingwidthblack = [125, 120, 160];
        stimpara.contrast = [1, 1, 1];
        stimpara.meanintensity = [0.5, 0.5, 0.5];
        stimpara.coneIsolating = false;
        stimpara.color = false;
        
    case {'directiongratingsequence_bw_2_stimuli_88angles','opto_directiongratingsequence_bw_2_stimuli_88angles'}
        stimpara.Nangles = [8,8];
        stimpara.numstimuli = 2;
        stimpara.cycles = [4, 4];
        stimpara.period = [100, 15];
        stimpara.duration = [400, 180];
        stimpara.regeneration = [300, 120];
        stimpara.preframestimuli = 180;
        stimpara.squarewave = [0, 1];
        stimpara.gratingwidthwhite = [125, 120];
        stimpara.gratingwidthblack = [125, 120];
        stimpara.contrast = [1, 1];
        stimpara.meanintensity = [0.5, 0.5, 0.5];
        stimpara.coneIsolating = false;
        stimpara.color = false;
        
    case {'directiongratingsequence_marmoset_6stimuli', 'directiongratingsequence_marmoset_6stim',...
            'directiongratingseq_marmoset_6stimuli','directiongratingsequence_bw_6_stimuli_marmoset'}
        stimpara.color = false;
        stimpara.numstimuli = 6;
        stimpara.Nangles = [8 8 8 8 8 8];
        stimpara.period = [100 100 30 30 15 15];
        stimpara.squarewave = [0 0 0 0 0 0];
        stimpara.cycles = [4 4 4 4 4 4];
        stimpara.duration = [300 300 300 300 300 300 ];
        stimpara.regeneration = [60 60 60 60 60 60];    %preframes Angles
        stimpara.preframestimuli = 0;
        stimpara.gratingwidthwhite = [120 240 120 240 120 240];
        stimpara.gratingwidthblack = [120 240 120 240 120 240];
        stimpara.contrast = [1 1 1 1 1 1];
        stimpara.meanintensity = [0.5 0.5 0.5 0.5 0.5 0.5];
        stimpara.coneIsolating = false;
        
    case {'directiongratingsequence_marmoset_6stimuli_2to5hz','directiongrating_marmoset_6stimuli_2to5hz'...
            'directiongratingsequence_marmoset_6stimuli_2to5hz1','directiongratingseq_marmoset_6stimuli_2to5hz',...
            'directiongratingsequence_marmoset_6stimuli_2to5hz_after','dssequence_2to5hz'}
        stimpara.color = false;
        stimpara.numstimuli = 6;
        stimpara.Nangles = [8 8 8 8 8 8];
        stimpara.period = [30 30 15 15 12 12];
        stimpara.squarewave = [0 0 0 0 0 0];
        stimpara.cycles = [4 4 4 4 4 4];
        stimpara.duration = [240 240 240 240 240 240];
        stimpara.regeneration = [60 60 60 60 60 60];    %preframes Angles
        stimpara.preframestimuli = 0;
        stimpara.gratingwidthwhite = [120 180 120 180 120 180];
        stimpara.gratingwidthblack = [120 180 120 180 120 180];
        stimpara.contrast = [1 1 1 1 1 1];
        stimpara.meanintensity = [0.5 0.5 0.5 0.5 0.5 0.5];
        stimpara.coneIsolating = false;
        
    case 'movingbars8dirs_bkg0_speed900'
        clearvars stimpara;
        stimpara.stimulus = 'movingbar';
        stimpara.height = 4005;
        stimpara.length = 600;
        stimpara.speed = 900;
        stimpara.preframes = 300;
        stimpara.angles = [0 45 90 135 180 225 270 315];
        stimpara.step = 0;
        stimpara.contrast = 1;
        stimpara.background = 0;
        stimpara.barintensity = 1;
        
    case {'plaidpatches8dir_150deg_period80_pxl_v1_', 'plaidpatches8dir_150deg_period80_pxl_v1'}
        clearvars stimpara;
        stimpara.stimulus = 'plaidpatches';
        stimpara.seed = -10000;
        stimpara.nblinks = 1;
        stimpara.stixelwidth = 1;
        stimpara.stixelheight = 1;
        stimpara.radius = 300;
        stimpara.velocity = 1;
        stimpara.period = 80;
        stimpara.ratio = 0.3;
        stimpara.angle = 150;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.5;
        stimpara.transparency = 0.5;
        stimpara.direction = 0;
        stimpara.ndir = 8;
        stimpara.repeats = 5;
        stimpara.preFrames = 100;
        stimpara.lmargin = 0;
        stimpara.rmargin = 0;
        stimpara.bmargin = 0;
        stimpara.tmargin = 0;
        
    case 'directiongratingsequence_2stimuli_ds_os_for75hz'
        stimpara.color = false;
        stimpara.numstimuli = 2;
        stimpara.Nangles = [8,8];
        stimpara.period = [125, 20];
        stimpara.squarewave = [0, 1];
        stimpara.cycles = [5, 5];
        stimpara.regeneration = [225, 150];
        stimpara.preframestimuli = 0;
        stimpara.duration = [500, 240];
        stimpara.gratingwidthwhite = [120, 112.5];
        stimpara.gratingwidthblack = [120, 112.5];
        stimpara.contrast = [1, 1];
        stimpara.meanintensity = [0.5, 0.5, 0.5];
        stimpara.coneIsolating = false;
        
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end
