

function [ newStimpara ] = stimCheckerFlickerPlusMovie (desc, stimpara)
%
%%% stimCheckerFlickerPlusMovie %%%
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

% defualt values
stimpara.stimulus = lower('checkerflickerplusmovie');
stimpara.stixelwidth = 8;
stimpara.stixelheight = 8;
stimpara.blackwhite = true;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.Nblinks = 2;
stimpara.Nblinksmovie = 2;
stimpara.seed = -1000;
stimpara.secondseed = -10000;
stimpara.runningframes = 1500;
stimpara.color = false;
if stimpara.color
    stimpara.independentcolors = true;
    stimpara.coneIsolating = true;
    stimpara.usered = false;
    stimpara.usegreen = true;
    stimpara.useblue = true;
    stimpara.redmean = 0.5;
    stimpara.greenmean = 0.5;
    stimpara.bluemean = 0.5;
    if stimpara.blackwhite, stimpara.redContrast = 0.75; else, stimpara.redContrast = 0.3; end
    if stimpara.blackwhite, stimpara.greenContrast = 0.75; else, stimpara.greenContrast = 0.3; end
    if stimpara.blackwhite, stimpara.blueContrast = 0.75; else, stimpara.blueContrast = 0.3; end
end
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;
stimpara.path = 'Y:\FromPeopleToPeople\Dimos\natural_stimuli\doves_mov1_subj5_60Hz\';
stimpara.width = 800;
stimpara.height = 600;


switch lower(desc)
    
    case {'checkerflickerplusmovie2x2bw2blinks_doves60hz','checkerflicker2x2bw2blinks_doves60hz1'}
        stimpara.nblinksmovie = 1;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        
    case 'checkerflickerplusmovie2x2bw3blinks_doves60hz'
        stimpara.nblinksmovie = 1;
        stimpara.runningframes = 1200;
        stimpara.frozenframes = 120;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 3;
        
    case 'checkerflickerplusmovie4x4bw2blinks_doves60hz'
        stimpara.nblinksmovie = 1;
        stimpara.runningframes = 1800;
        stimpara.frozenframes = 180;
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
        
    case {'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed1'}
        stimpara.nblinksmovie = 1;
        stimpara.runningframes = 1200;
        stimpara.frozenframes = 120;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 3;
        stimpara.seed = -1000;
        
    case {'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed2'}
        stimpara.nblinksmovie = 1;
        stimpara.runningframes = 1200;
        stimpara.frozenframes = 120;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 3;
        stimpara.seed = -1001;
        
    case {'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed3'}
        stimpara.nblinksmovie = 1;
        stimpara.runningframes = 1200;
        stimpara.frozenframes = 120;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.nblinks = 3;
        stimpara.seed = -1002;
        
end

if ~stimpara.color, stimpara = rmfield(stimpara,'color'); end


newStimpara = stimpara;

end
