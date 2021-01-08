

function [ newStimpara ] = stimTempFreqSensitivity (desc, stimpara)
%
%%% stimTempFreqSensitivity %%%
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
    
    case 'tempfrequsensitivity_contrast0p8'
        stimpara.contrast = 0.8;
        
    case 'tempfrequsensitivity_contrast0p1'
        stimpara.contrast = 1;
        
    case {'tempfrequsensitivity120preframes60levels8n', 'tempfrequsensitivity120preframes60levels8'...
            'tempfrequesensivity120preframes60levels8','tempfrequencysensitivity120preframes60levels8n',...
            'tempfrequsensitivity120preframes60level8n','tempfreqsensitivity120preframes60levels8n',...
            'tempfreqsensitivity120preframes60levels8'}
        stimpara.nframes = 120;
        stimpara.preframes = 60;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 0.5;
        stimpara.nfrequencies = 12;
        %stimpara.freqs = [1,2,3,5,7,9,12,15];
        stimpara.frequencylist = [1,2,3,5,7,9,12,15];
        
    case 'tempfrequsensitivity_contrast0p8_lap4'
        stimpara.contrast = 0.8;
        stimpara.pharmacology = true;
        stimpara.drugName = 'LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'tempfrequsensitivity_contrast0p8_cnqx'
        stimpara.contrast = 0.8;
        stimpara.pharmacology = true;
        stimpara.drugName = 'CNQX and LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'tempfrequsensitivity_contrast0p8_washout'
        stimpara.contrast = 0.8;
        stimpara.pharmacology = true;
        stimpara.drugName = 'washout';
        stimpara.drugDose = 'Ames solution';
end

newStimpara = stimpara;

end


