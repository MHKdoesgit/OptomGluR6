

function [ newStimpara ] = stimONOFFsteps (desc, stimpara)
%
%%% stimONOFFsteps %%%
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
    
    case {'onoffsteps30'
            'onoffsteps_contrast1_30frames'
            'onoffsteps30_50'
            'onoffstep30_50'}
        stimpara.nframes = 30;
        
    case 'onoffsteps30_contrast0.3'
        stimpara.nframes = 30;
        stimpara.contrast = 0.3;
        
    case 'onoffsteps30_preframes90_contrast0.4'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.contrast = 0.4;
        
    case {'onoffsteps30preframes90contrast1p', 'onoffsteps30preframes90contrasts1',...
            'onoffsteps30_preframes90_contrast1p','onoffsteps30preframes90contrast1p1',...
            'onoffsteps30preframes90contrast1','onoffsteps3w0preframes90contrast1p',...
            'onoffsteps30preframes90contrast1_control','onoffsteps30preframes90contrast1_am251_',...
            'onoffstep30preframes90contrast1_after','onoffsteps30preframes90contrast1_win55212',...
            'onoffsteps30preframes90contrast1_after','onoffsteps30preframes90contrast1_before1h'}
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.contrast = 1;
        
    case 'onoffsteps_60blinks_60preframes'
        stimpara.nframes = 60;
        stimpara.preframes = 60;
        stimpara.contrast = 1;
    
    case 'onoffsteps_60blinks'
        stimpara.nframes = 60;
        stimpara.preframes = 0;
        stimpara.contrast = 1;
        
    case 'onoffsteps30preframes60contrast1p'
        stimpara.nframes = 30;
        stimpara.preframes = 60;
        stimpara.contrast = 1;
        
    case 'onoffstep30preframes90contrast1_bright79'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.contrast = 1;
        stimpara.OLEDbrightness = 79;
        
    case   'onoffsteps30preframes90contrast1_bright91'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.contrast = 1;
        stimpara.OLEDbrightness = 91;
        
    case {'onoffsteps10_preframes30_contrast0.5'...
            'onoffsteps10_preframes30_0.5contrast'}
        stimpara.nframes = 10;
        stimpara.preframes = 30;
        stimpara.contrast = 0.5;
        
    case {'onoffstep30_50contrast'
            'onoffsteps30_50contrast'}
        stimpara.nframes = 30;
        stimpara.preframes = 0;
        stimpara.contrast = 1;
        
    case {'onoffsteps30_preframes90_contrast1','onoffsteps_100contrast_30frames_90pfr'}
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.contrast = 1;
        
    case {'onoffsteps_contrast1_30frames120preframes', 'onoffsteps30_preframes120_contrast1'}
        stimpara.nframes = 30;
        stimpara.preframes = 120;
        stimpara.contrast = 1;
        
    case {'onoffsteps_100contrast_60frames_90pfr','opto_onoffsteps_100contrast_60frames_90pfr'}
        stimpara.nframes = 60;
        stimpara.preframes = 90;
        stimpara.contrast = 1;
        
    case {'onoffsteps_100contrast_12frames_288pfr'}
        stimpara.nframes = 12;
        stimpara.preframes = 288;
        stimpara.contrast = 1;
        
    case 'onoffsteps30_preframes120_contrast0.6'
        stimpara.nframes = 30;
        stimpara.preframes = 120;
        stimpara.contrast = 0.6;
        
    case 'onoffsteps30_preframes90_contrast1_lap4'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.pharmacology = true;
        stimpara.drugName = 'LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'onoffsteps30_preframes90_contrast1_cnqx'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.pharmacology = true;
        stimpara.drugName = 'CNQX and LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'onoffsteps30_preframes90_contrast1_washout'
        stimpara.nframes = 30;
        stimpara.preframes = 90;
        stimpara.pharmacology = true;
        stimpara.drugName = 'washout';
        stimpara.drugDose = 'Ames solution';
        
    case {'onoffsteps_3600_for_iprgcs', 'onoffsteps3600_for_iprgc','opto_onoffsteps3600_for_iprgc'}
        stimpara.nframes = 3600;
        stimpara.preframes = 0;
        stimpara.contrast = 1;
        
    case {'onoffsteps1800_for_iprgc'}
        stimpara.nframes = 1800;
        stimpara.preframes = 0;
        stimpara.contrast = 1;
        
    case 'onoffsteps1200_preframes600_contrast0p6_iprgcs'
        stimpara.nframes = 1200;
        stimpara.preframes = 600;
        stimpara.contrast = 0.6;
        
    case {'onoffsteps75_preframe150_contrast1_for75hz','onoffsteps75_preframes150_contrast1_for75hz'}
        stimpara.nframes = 75;
        stimpara.preframes = 150;
        stimpara.contrast = 1;
end

newStimpara = stimpara;
end