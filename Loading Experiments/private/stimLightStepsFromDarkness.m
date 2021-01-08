

function [ newStimpara ] = stimLightStepsFromDarkness (desc, stimpara)
%
%%% stimLightStepsFromDarkness %%%
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

% default values
stimpara.stimulus = 'lightstepsfromdarkness';
stimpara.nframes = 60;
stimpara.preframes = 240;
stimpara.nsteps = 10;
stimpara.nrepeats = 2;
stimpara.nintensities = 8;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;


switch lower(desc)
    
    case {'lightstepsfromdarkness60pfr240with10steps'
            'lightstepsfromdarkness60pfr240with10steps_washout'}
        
        if strcmp(desc,'lightstepsfromdarkness60pfr240with10steps_washout')
            stimpara.pharmacology = true;
            stimpara.drugName = 'washout';
            stimpara.drugDose = 'Ames solution';
        end
        
    case {'lightstepsfromdarkness_60stim_240pfr_20steps_4repeats'
            'opto_lightstepsfromdarkness_60stim_240pfr_20steps_4rep'}
        stimpara.nrepeats = 4;
        
    case {'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats'
            'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats_40'
            'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats_0'}
        stimpara.nrepeats = 3;
        
    case {'lightstepsfromdarkness_75stim_300pfr_10steps_3repeats'}
        stimpara.nframes = 75;
        stimpara.preframes = 300;
        stimpara.nsteps = 10;
        stimpara.nrepeats = 3;
end

newStimpara = stimpara;


end

