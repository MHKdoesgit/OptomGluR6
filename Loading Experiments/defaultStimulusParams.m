function [ stimpara ] = defaultStimulusParams( desc, stimpara, screenWidth, screenHeight )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch lower(desc)
    case 'barcodestimulus'
        stimpara = defaultBarcodeStimulus(stimpara);
    case 'checkerflicker'
        stimpara = defaultCheckerFlicker(stimpara);
    case 'checkerflickerplusmovie'
        stimpara = defaultCheckerFlickerPlusMovie(stimpara);
    case 'chirpstimulus'
        stimpara = defaultChirpStimulus(stimpara);
    case 'contrastsensitivitywithdriftinggratings'
        stimpara = defaultContrastSensitivityWithDriftingGratings(stimpara);
    case 'contraststepladder'
        stimpara=defaultContrastStepLadder(stimpara);
    case 'directiongrating'
         stimpara = defaultDirectionGrating(stimpara);
    case 'directiongratingsequence'
        stimpara = defaultDirectionGratingSequence(stimpara);
    case 'frozennoise'
        stimpara = defaultFrozenNoise(stimpara);
    case 'fixationmovie'
        stimpara = defaultFixationMovie(stimpara);   
    case 'gratingflashes'
        stimpara = defaultGratingFlashes(stimpara);   
    case 'imagesequence'
        stimpara = defaultImageSequence(stimpara);
    case 'locallysparsenoise'
        stimpara = defaultLocallySparseNoise(stimpara);
    case 'movingbars'
        stimpara = defaultMovingBars(stimpara);
    case 'onoffgrating'
        %no defaults are specified by the program
    case 'onoffsteps'
        stimpara=defaultOnOffSteps(stimpara);
    case 'rotatingstripes'
        stimpara=defaultRotatingStripes(stimpara,  screenWidth, screenHeight);
    case 'reversinggratingswithvaryingspatialperiod'
        stimpara=defaultReversingGratingsWithVaryingSpatialPeriod(stimpara);
    case 'reversinggrating'
        stimpara=defaultReversingGrating(stimpara);
    case 'locallysparsesubunitflash'
        stimpara=defaultLocallySparseSubunitFlash(stimpara);
    case 'sparsesubunitflash'
        stimpara=defaultSparseSubunitFlash(stimpara);
    case 'saccadegrating'	
        stimpara=defaultSaccadeGrating(stimpara);
    case 'subunitflash'
        stimpara=defaultSubunitFlash(stimpara);
    case 'tempfrequsensitivity'
        stimpara=defaultTempFrequSensitivity(stimpara);
    case 'triggereddriftinggratings'
        stimpara=defaultTriggeredDriftingGratings(stimpara);
        
end

end

function stimpara = defaultBarcodeStimulus(stimpara)

stimpara.preframes     = 60;
stimpara.speed         = 2;
stimpara.contrast      = 1;
stimpara.maxperiod     = 1280;
stimpara.minperiod     = 16;
stimpara.barcodelength = 1280;
stimpara.usevertical   = true;
stimpara.meanintensity = 0.5;
stimpara.lmargin       = 0;
stimpara.rmargin       = 0;
stimpara.bmargin       = 0;
stimpara.tmargin       = 0;
stimpara.seed          = -10000;

end

function stimpara= defaultCheckerFlicker(stimpara)
stimpara.stixelwidth = 10;
stimpara.stixelheight = 10;
stimpara.blackwhite = true;
stimpara.meanintensity = 0.5;

stimpara.contrast = 1;
stimpara.Nblinks = 1;
stimpara.pulseRate = 2;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

stimpara.seed = -10000;
stimpara.secondSeedFlag=false;

%added to give equivalence to frozen stimuli
%stimpara.RunningFrames=0;
%stimpara.FrozenFrames=0;

stimpara.color = false;
stimpara.useRed = false;
stimpara.useGreen = true;
stimpara.useBlue = true;    
stimpara.redMeanIntensity = 0.5;
stimpara.greenMeanIntensity = 0.5;
stimpara.blueMeanIntensity = 0.5;
end

function stimpara = defaultContrastSensitivityWithDriftingGratings(stimpara)
stimpara.Nrepeats       = 1;
stimpara.Nframes        = 180;
stimpara.preframes      = 60;
stimpara.postframes     = 60;
stimpara.meanintensity  = 0.5;
stimpara.temporalPeriod = 30;
end

function stimpara=defaultContrastStepLadder(stimpara)
stimpara.Nframes=12;
stimpara.preframes=108;
stimpara.meanintensity=0.5;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
stimpara.contrastlist=[0 0.01 0.02 0.03 0.05 0.08 0.12 0.18 0.24];
end



function stimpara=defaultFrozenNoise(stimpara)
stimpara.stixelwidth = 8;
stimpara.stixelheight = 8;
stimpara.blackwhite = true;
stimpara.meanintensity = 0.5;

stimpara.contrast = 1;
stimpara.Nblinks = 2;
stimpara.pulseRate = 2;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

stimpara.seed = -1000;
stimpara.secondseed = -10000;

stimpara.radius=60;
stimpara.RunningFrames=1000;
stimpara.FrozenFrames=1000;
stimpara.useMask=0;
stimpara.centerX=400;
stimpara.centerY=300;
end

function stimpara=defaultCheckerFlickerPlusMovie(stimpara)
stimpara.stixelwidth = 8;
stimpara.stixelheight = 8;
stimpara.blackwhite = true;
stimpara.meanintensity = 0.5;

stimpara.contrast = 1;
stimpara.Nblinks = 2;
stimpara.pulseRate = 2;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

stimpara.seed = -1000;
stimpara.secondseed = -10000;

stimpara.radius=60;
stimpara.RunningFrames=1500;
stimpara.FrozenFrames=0;

stimpara.path='Y:\FromPeopleToPeople\Dimos\natural_stimuli\';
stimpara.Nblinksmovie=2;
stimpara.width=800;
stimpara.height=600;
end

function stimpara=defaultChirpStimulus(stimpara)
stimpara.color=0;
stimpara.Nrepeats=10;
stimpara.meanintensity=0.5;
stimpara.maxcontrast=0.3;
stimpara.preframes=120;

stimpara.onoffstep.duration=3*60;
stimpara.onoffstep.preframes=2*2*60;
stimpara.onoffstep.contrast=1;
stimpara.onoffstep.preframesintensity=0;

stimpara.freqsweep.preframes=2*60;
stimpara.freqsweep.duration=8*60;
stimpara.freqsweep.lo_freq=0;
stimpara.freqsweep.hi_freq=8;
stimpara.freqsweep.contrast=1;

stimpara.ctrsweep.preframes=2*60;
stimpara.ctrsweep.duration=8*60;
stimpara.ctrsweep.freq=2;
stimpara.ctrsweep.contrast=1;

stimpara.lmargin=0;
stimpara.rmargin=0;
stimpara.bmargin=0;
stimpara.tmargin=0;

stimpara.useMask=0;
stimpara.centerX=400;
stimpara.centerY=300;
stimpara.radius=60;
end

function stimpara= defaultDirectionGrating(stimpara)
stimpara.seed = -10000;
stimpara.period = 80;
stimpara.cycles=3;
stimpara.sequence=false;
stimpara.regeneration = 100;
stimpara.brightness = 1;
stimpara.Nangles = 8;
stimpara.duration = 400;
stimpara.gratingwidth = 350;
end

function stimpara= defaultDirectionGratingSequence(stimpara)
stimpara.color = true;
stimpara.nangles = 8;
stimpara.period = [100 60 60 150 70 15];
stimpara.squareWave=[0 1 1 0 0 1];
stimpara.cycles=[5 5 5 5 5];
stimpara.preframeAngles = [ 300, 120, 120, 300, 120, 120 ];
stimpara.preframeStimuli = 180;
stimpara.gratingwidthwhite = [ 117.5  300 425 120 141 112.5 ];
stimpara.gratingwidthblack = [ 117.5  900 425 120 141 112.5 ];
stimpara.contrasts = [1 1 1 1 1 1];
stimpara.meanintensity = [0.5 0.5 0.5];
stimpara.duration =  [200 240 240 600 280 180];
end

function stimpara = defaultFixationMovie(stimpara)

stimpara.meanintensity    = 0.5;
stimpara.FrozenFixations  = 0;
stimpara.RunningFixations = 1500;
stimpara.lmargin          = 0;
stimpara.rmargin          = 0;
stimpara.bmargin          = 0;
stimpara.tmargin          = 0;

end

function stimpara = defaultGratingFlashes(stimpara)

stimpara.meanintensity = 0.5;
stimpara.stripewidths  = [1, 2,  4,  8, 16, 32, 64, 128];
stimpara.Norientations = 10;
stimpara.Nphases       = 4;
stimpara.lmargin       = 0;
stimpara.rmargin       = 0;
stimpara.bmargin       = 0;
stimpara.tmargin       = 0;
stimpara.seed          = -10000;
stimpara.stimduration  = 45;
stimpara.flashstart    = 10;
stimpara.flashstop     = 19;
stimpara.nrepeats      = 10;

end

function stimpara=defaultImageSequence(stimpara)
stimpara.path='Y:\FromPeopleToPeople\Dimos\natural_stimuli';
stimpara.Nx=800;
stimpara.Ny=600;
stimpara.width=1;
stimpara.height=1;
stimpara.trialduration=60;
stimpara.israndom=1;
stimpara.nrepeats=10;
stimpara.backgroundintensity=0.5;
stimpara.seed = -10000;
end

function stimpara=defaultLocallySparseNoise(stimpara)
stimpara.stixelwidth = 10;
stimpara.stixelheight = 10;
stimpara.gapwidth = 2;
stimpara.gapheight = 2;
stimpara.meanintensity=0.5;
stimpara.contrast=1;
stimpara.stimduration=2;
stimpara.seedlocation= -1000;
stimpara.seedcontrast = -10000;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
end

function stimpara = defaultMovingBars(stimpara)

stimpara.height         = 40;
stimpara.length         = 120;
stimpara.horizontalgap  = 4 * 120;
stimpara.verticalgap    = 10;
stimpara.Nbars          = 5;
stimpara.maskdiameter   = 800;
stimpara.speed          = 2;
stimpara.preframes      = 0;
stimpara.Nangles        = 8;
stimpara.usedark        = true;
stimpara.contrast       = 1;
stimpara.background     = 0.5;
        
end

function stimpara=defaultSparseSubunitFlash(stimpara)
stimpara.stixelside=14;
stimpara.nstixels=2;
stimpara.npositionsx=8;
stimpara.npositionsy=8;
stimpara.ndist=2;
stimpara.showtwo=0;
stimpara.nangles=8;
stimpara.ncontrasts=5;
stimpara.mincontrast=0.2;
stimpara.maxcontrast=1;
stimpara.stimduration=12;
stimpara.nrepeats=4;
stimpara.meanintensity=0.5;
stimpara.seed=-10000;
end


function stimpara=  defaultOnOffSteps(stimpara)
stimpara.Nframes = 60;
stimpara.preframes = 0;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.tmargin = 0;
stimpara.bmargin = 0;
end

function stimpara=defaultReversingGratingsWithVaryingSpatialPeriod(stimpara)
stimpara.Nframes=60;
stimpara.Nreversals=20;
stimpara.preframes=120;
stimpara.Nrepeats=2;
stimpara.contrast=1;
stimpara.meanintensity=0.5;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
stimpara.grayframes=0;
stimpara.stripewidths=[2 4 8 16 32];
stimpara.Nphases=[2 2 2 2 2];
end

function stimpara=defaultReversingGrating(stimpara)
stimpara.stripewidth= 10;
stimpara.Nframes = 30;
stimpara.grayframes= 0;
stimpara.contrast = 1;
stimpara.meanintensity = 0.5;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
end

function stimpara=defaultRotatingStripes(stimpara,  screenWidth, screenHeight)
stimpara.stixelwidth=10;
stimpara.blackwhite=true;
stimpara.meanintensity=0.5;
stimpara.contrast = 1;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

stimpara.color = false;

stimpara.stimduration=600;
stimpara.preframes=120;
stimpara.nangles=4;
stimpara.seed=-10000;

stimpara.xoffset=0;
stimpara.yoffset=0;

stimpara.Nblinks=1;
stimpara.maskdiameter=min([screenWidth screenHeight ]);
end

function stimpara=defaultSaccadeGrating(stimpara)

stimpara.fixationframes=80;
stimpara.saccadeframes=10;
stimpara.meanintensity=0.5;
stimpara.highintensity = 0.8;

stimpara.barwidth = 40;
stimpara.cosineprofile = false;
stimpara.averageshift = 2;
stimpara.centralgrey = false;
stimpara.contrastchanges = false;
stimpara.highresolution = false;

stimpara.centralwidth=1100;
stimpara.centralheight=1100;
stimpara.centreOffsetX=0;
stimpara.centreOffsetY=0;
stimpara.seed = -10000;
end


function stimpara = defaultSubunitFlash(stimpara)

stimpara.stixelwidth    = 10;
stimpara.stixelheight   = 10;
stimpara.meanintensity  = 0.5;
stimpara.nangles        = 8;
stimpara.ncontrasts     = 8;
stimpara.mincontrast    = 0;
stimpara.maxcontrast    = 1;

stimpara.lmargin        = 0;
stimpara.rmargin        = 0;
stimpara.bmargin        = 0;
stimpara.tmargin        = 0;

stimpara.seed           = -10000;
stimpara.stimduration   = 120;
stimpara.flashstart     = 30;
stimpara.flashstop      = 42;
stimpara.nrepeats       = 10;
stimpara.squaresampling = false;

end

function stimpara=defaultTempFrequSensitivity(stimpara)
stimpara.Nframes=300;
stimpara.preframes=120;
stimpara.meanintensity=0.5;
stimpara.contrast=1;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
stimpara.frequencylist=1:12;
end

function stimpara=defaultTriggeredDriftingGratings(stimpara)
stimpara.Nrepeats = 6;
stimpara.Nframes = 180;
stimpara.preframes = 60;
stimpara.postframes = 60;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.radius=10;
stimpara.temporalPeriod=20;
stimpara.centerX=0;
stimpara.centerY=0;
stimpara.untriggeredMode=false;
stimpara.spatialperiod=[1730 173 109.5 68.9 43.5 27.4 17.3 10.95 6.89 4.35];
end


function stimpara=defaultLocallySparseSubunitFlash(stimpara)
stimpara.stixelwidth = 10;
stimpara.stixelheight = 10;
stimpara.nstixelsx = 2;
stimpara.nstixelsy = 2;
stimpara.gapwidth = 2;
stimpara.gapheight = 2;

stimpara.nangles=8;
stimpara.ncontrasts=5;
stimpara.mincontrast=0.2;
stimpara.maxcontrast=1;
stimpara.stimduration=12;


stimpara.meanintensity=0.5;
stimpara.stimduration=12;

stimpara.seedlocation= -1000;
stimpara.seedcombination = -10000;

stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
end


