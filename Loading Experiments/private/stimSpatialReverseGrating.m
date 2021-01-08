

function [ newStimpara ] = stimSpatialReverseGrating (desc, stimpara)
%
%%% stimSpatialReverseGrating %%%
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
stimpara.stimulus = lower('reverse_grating_with_varing_spatial_period');
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)
    
    case 'reversinggratingwithvaryingspatialperiod'
        stimpara.nframes = 60;
        stimpara.preframes = 120;
        stimpara.nreversals = 20;
        stimpara.nrepeats = 2;
        
    case 'reversinggratingwvsp_30reversals_6widths'
        stimpara.nframes = 60;
        stimpara.preframes = 120;
        stimpara.nreversals = 30;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [2 4 8 14 28 1600];
        stimpara.nphases = [1 2 2 4 4 1];
        
    case 'reversinggratingwvsp12preframes60'
        stimpara.nframes = 12;
        stimpara.preframes = 60;
        stimpara.nreversals = 25;
        stimpara.grayframes = 0;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [1 2 4 8 16 32 64 1000];
        stimpara.nphases = [1 1 2 2 4 4 8 1];
        
    case {'reversinggratingwvsp_30reversals_24816321000_6widths','opto_reversinggratingwvsp_30reversals_24816321000_6widths',...
            'opto_reversinggratingwvsp_30rev_24816321000_6widths'}
        stimpara.nframes = 60;
        stimpara.preframes = 120;
        stimpara.nreversals = 30;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [2 4 8 16 32 1000];
        stimpara.nphases = [1 2 2 4 4 1];
        
    case {'reversinggratingwvsp12preframes60n','reversinggratingwvp12preframes60'...
            'reversingratingwvsp12preframes60', 'reversinggratingwvsp12preframes60f',...
            'reversingratingwvsp12preframes60f','reversinggratingwvsp_12stim_25reversals_8widths',...
            'reversinggratingwvsp12preframes'} % marmoset version
        stimpara.nframes = 12;
        stimpara.preframes = 60;
        stimpara.nreversals = 25;
        stimpara.grayframes = 0;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [1 2 4 8 16 32 64 800];
        stimpara.nphases = [1 1 2 2 4 4 8 1];
        
    case {'reversinggratingwvsp30preframes'}
        stimpara.nframes = 30;
        stimpara.preframes = 60;
        stimpara.nreversals = 25;
        stimpara.grayframes = 0;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [1 2 4 8 16 32 64 800];
        stimpara.nphases = [1 1 2 2 4 4 8 1];
        
    case {'reversinggratingwvsp45preframes90_for75hz_40br'
            'reversinggratingwvsp45preframes90_for75hz_0br'}
        stimpara.nframes = 45;
        stimpara.preframes = 90;
        stimpara.nreversals = 25;
        stimpara.grayframes = 0;
        stimpara.nrepeats = 2;
        stimpara.stripewidth = [1 2 4 8 16 32 64 800];
        stimpara.nphases = [1 1 2 2 4 4 8 1];
        
end

newStimpara = stimpara;

end

