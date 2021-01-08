

function [ newStimpara ] = stimImageSequence (desc, stimpara)
%
%%% stimImageSequence %%%
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
stimpara.stimulus = lower('imagesequence');
stimpara.trialduration = 90;
stimpara.flashstart = 30;
stimpara.flashstop = 42;
stimpara.nrepeats = 10;
stimpara.israndom = false;
stimpara.seed = -10000;
stimpara.backgroundIntensity = 0.5;
stimpara.width = 1;
stimpara.height = 1;
stimpara.imagepath = [];


switch lower(desc)
    
    case {'imgseq_512x512_300im_10r_mixed'
            'imgseq512x512_300im_10r_mixed'}
        
        stimpara.path = 'Y:\FromPeopleToPeople\Dimos\natural_stimuli\imgmix512x512_300\';
        stimpara.Nx = 512;
        stimpara.Ny = 512;
        stimpara.trialduration = 60;
        stimpara.flashstart = 8;
        stimpara.flashstop = 20;
        stimpara.israndom = true;
        stimpara.nrepeats = 10;
        
    case {'imgseq512x512_40im_10r_scales5_15to240um'
            'imgseq512x512_40im_10r_scale5_15to240p'
            'imgseq512x512_40im_10r_scale5_15_to_240'}
        
        stimpara.path = 'Y:\FromPeopleToPeople\Dimos\natural_stimuli\imgmix512x512_40_blur_scales5_15to240\';
        stimpara.Nx = 512;
        stimpara.Ny = 512;
        stimpara.trialduration = 60;
        stimpara.flashstart = 8;
        stimpara.flashstop = 20;
        stimpara.israndom = true;
        stimpara.nrepeats = 10;
        
    case {'imgseq512x512_40im_10r_blur_scales5_30to180'}
        stimpara.path = 'Y:\FromPeopleToPeople\Dimos\natural_stimuli\imgmix512x512_40_blur_scales5_30to180\';
        stimpara.Nx = 512;
        stimpara.Ny = 512;
        stimpara.trialduration = 60;
        stimpara.flashstart = 8;
        stimpara.flashstop = 20;
        stimpara.israndom = true;
        stimpara.nrepeats = 10;
        
    case  {'imgseq_movies800x600_doves30hz_blur_scales5_15_to240_20r','imgseq_movie800x600_doves30hz_blur_scales5'}
        
        stimpara.path = 'Y:\FromPeopleToPeople\Dimos\natural_stimuli\doves_mov1_subj5_30Hz_blur_scales5_15to240\';
        stimpara.Nx = 800;
        stimpara.Ny = 600;
        stimpara.trialduration = 2;
        stimpara.flashstart = 0;
        stimpara.flashstop = 2;
        stimpara.israndom = false;
        stimpara.nrepeats = 20;
end

newStimpara = stimpara;

end
