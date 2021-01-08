

function [ newStimpara ] = loadDefault(stimName, stimpara)
%
%%% loadDefault %%%
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


if nargin < 2
    error('not enough arguments (need 2)');
end

stimName = lower(stimName);
stimpara.stimulus = stimName;

switch stimName
    
    case 'checkerflicker'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.blackwhite = true;
        stimpara.contrast = 1;
        stimpara.color = false;
        stimpara.meanintensity = 0.5;
        stimpara.seed = -1000;
        stimpara.secondseed = -10000;
        stimpara.thirdseed = -2000;
        stimpara.nblinks = 2;
        stimpara.coneisolating = false;
        
    case 'fullfieldflicker'
        stimpara.stixelwidth = 1600;
        stimpara.stixelheight = 1200;
        stimpara.blackwhite = false;
        stimpara.color = false;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.3;
        stimpara.seed = -1000;
        stimpara.nblinks = 2;
        stimpara.coneisolating = false;
        
    case 'stripsflicker'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 600;
        stimpara.blackwhite = true;
        stimpara.contrast = 1;
        stimpara.color = false;
        stimpara.meanintensity = 0.5;
        stimpara.seed = -10000;
        stimpara.nblinks = 2;
        
    case 'pinknoise'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.meanintensity = 0.5;
        stimpara.contrastwhite = 0.3;
        stimpara.contrastpink = 0.3;
        stimpara.fracwhite = 0.0;
        stimpara.beta = 2.0;
        stimpara.temporalbeta = 1.0;
        stimpara.nodes = 300;
        stimpara.seed = -10000;
        stimpara.nblinks = 2;
        
    case 'onoffsteps'
        stimpara.nframes = 30;
        stimpara.preframes = 0;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        
    case 'onoffgrating'
        stimpara.squarewave = false;
        stimpara.period = 60;
        stimpara.speed = 30;
        stimpara.phaseduration = 30;
        stimpara.angle = 0;
        
    case 'directiongrating'
        stimpara.seed = -10000;
        stimpara.regeneration = 0;
        stimpara.Nangles = 8;
        stimpara.meanintensity = 0.5;
        stimpara.period = 80;
        stimpara.duration = 400;
        stimpara.cycles = 10;
        stimpara.gratingwidth = 300;
        
    case 'contrastadaptation'
        stimpara.stixelwidth = 1600;
        stimpara.stixelheight = 1200;
        stimpara.seed = -1000;
        stimpara.blackwhite = false;
        stimpara.Nblinks = 2;
        stimpara.switchFrames = 1200;
        stimpara.highcontrast = 0.32;
        stimpara.lowcontrast = 0.08;
        
    case 'spatialcontrastadaptation'
        stimpara.stixelwidth = 8;
        stimpara.stixelheight = 8;
        stimpara.graywidth = 8;
        stimpara.grayheight = 8;
        stimpara.seed = -1000;
        stimpara.blackwhite = true;
        stimpara.Nblinks = 2;
        stimpara.switchFrames = 1200;
        stimpara.highcontrast = 1.0;
        stimpara.lowcontrast = 0.2;
        stimpara.twofixed = true;
        
    case 'rgbsteps'
        stimpara.nframes = 30;
        stimpara.preframes = 120;
        stimpara.meanintensity = 0.5; %for all color
        stimpara.contrast = 1; %for all contrast
        stimpara.color = true;
        stimpara.opponency = false;
        stimpara.isoResponse = false;
        stimpara.seedflag = true;
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redcontrast = 0;
        stimpara.greencontrast = 0.7;
        stimpara.bluecontrast = 0.7;
        
    case 'colorintensitysteps'
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmax = false;
        stimpara.bluemax = false;
        stimpara.greenmax = false;
        stimpara.redvalue = 1/256;
        stimpara.greenvalue = 1/256;
        stimpara.bluevalue = 1/256;
        stimpara.offsetvalue = 0;
        stimpara.framesonscreen = 60;
        stimpara.backgroundonscreen = 120;
        stimpara.skipvalue = 16;
        stimpara.nblinks = 2;
        stimpara.coneisolating = true;
        stimpara.backgroundIntensity = 128;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.maxcontrsat = 1;
        
    case 'rotatingstrips'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 600;
        stimpara.blackwhite = true;
        stimpara.contrast = 1;
        stimpara.color = false;
        stimpara.meanintensity = 0.5;
        stimpara.seed = -1000;
        stimpara.nblinks = 2;
        stimpara.stimduration = 10*60;
        stimpara.nAngles = 4;
        stimpara.coneisolating = false;
        stimpara.preframes = 120;
        
    case 'contraststepladder'
        stimpara.nframes = 12;
        stimpara.preframes = 180;
        stimpara.meanintensity = 0.5;
        stimpara.color = false;
        stimpara.colormode = 'white';
        stimpara.coneisolating = false;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'silentexchange'
        stimpara.nframes = 30;
        stimpara.preframes = 120;
        stimpara.meanintensity = 0.5; %for black & white
        stimpara.contrast = 1; %for black & white
        stimpara.color = true;
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redcontrast = 0;
        
    case 'coneisolationtest'
        stimpara.nframes = 30;
        stimpara.preframes = 120;
        stimpara.meanintensity = 0.5; %for black & white
        stimpara.contrast = 1; %for black & white
        stimpara.color = true;
        stimpara.coneisolating = false;
        stimpara.usered = false;
        stimpara.usegreen =true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redcontrast = 0;
        stimpara.bluecontrast = 0.7;
        
    case 'colorintegration'
        stimpara.stimduration = 30;
        stimpara.preframes = 120;
        stimpara.contrastdiff = 0.02;
        stimpara.mincontrast = -0.2; % -20% contrast
        stimpara.maxcontrast = 0.2; % 20% contrast
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
        
    case 'colorintegrationgrating'
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        stimpara.squarewave = false;
        stimpara.contrastdiff = 0.02;
        stimpara.mincontrast = -0.2; % -20% contrast
        stimpara.maxcontrast = 0.2; % 20% contrast
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorisoresponse'
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.numrepeats  = 3;
        stimpara.numangles = 36;
        stimpara.numcontrasts = 18;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.72; % to 72% contrast
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorisoresponsegrating'
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        stimpara.squarewave = false;
        stimpara.numrepeats  = 1;
        stimpara.numangles = 24;
        stimpara.numcontrasts = 9;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.72; % to 72% contrast
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'colorexchange'
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.contrastdiff = 0.20;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4; % 40% contrast
        stimpara.fixedcontrast = 0.20; % 10% fixed contrast
        stimpara.seed = -1000;
        stimpara.coneisolating = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        stimpara.redcontrast = 0;
        stimpara.greencontrast = 0.5;
        stimpara.bluecontrast = 0.5;
        
    case 'chirp'
        stimpara.Nrepeats = 100;
        stimpara.meanintensity = 0.5;
        stimpara.onoffstep.duration = 30;
        stimpara.onoffstep.preframes = 120;
        stimpara.onoffstep.contrast = 1;
        stimpara.freqsweep.preframes = 120;
        stimpara.freqsweep.duration = 300;
        stimpara.freqsweep.low_freq = 0.5;
        stimpara.freqsweep.high_freq = 8.0;
        stimpara.freqsweep.contrast = 1;
        stimpara.contrastsweep.preframes = 120;
        stimpara.contrastsweep.duration = 300;
        stimpara.contrastsweep.freq = 5.0;
        stimpara.contrastsweep.contrast = 1;
        
    case 'frozennoise'
        stimpara.stixelwidth = 1600;
        stimpara.stixelheight = 1200;
        stimpara.RunningFrames = 1500;
        stimpara.FrozenFrames = 300;
        stimpara.blackwhite = false;
        stimpara.color = false;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.3;
        stimpara.seedrunningnoise = -1000;
        stimpara.seedfrozennoise = -10000;
        stimpara.nblinks = 2;
        
    case 'saccadegrating'
        stimpara.fixationframes = 80;
        stimpara.saccadeframes = 10;
        stimpara.meanintensity = 0.5;
        stimpara.highintensity = 0.8;
        stimpara.lowintensity = 0.2;
        stimpara.barwidth = 40;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
        stimpara.centralgrey = false;
        stimpara.surroundgrey = false;
        stimpara.contrastchanges = false;
        stimpara.highresolution = false;
        stimpara.centralwidth = 1100;
        stimpara.centralheight = 1100;
        stimpara.centreOffsetX = 0;
        stimpara.centreOffsetY = 0;
        
    case 'tempfreqsensitivity'
        stimpara.nframes = 300;
        stimpara.preframes = 120;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara.nfrequencies = 12;
        stimpara.frequencylist = [1,2,3,4,5,6,7,8,9,10,11,12];
        
end

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

newStimpara = stimpara;

end

