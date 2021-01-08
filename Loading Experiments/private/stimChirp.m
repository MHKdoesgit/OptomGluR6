

function [ newStimpara ] = stimChirp (desc, stimpara)
%
%%% stimChirp %%%
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
    
    case {'chirpnf30pf120frq300lf05hf8ctr300frq5bluegreen'}
        stimpara.color = true;
        stimpara.onoffstep.redcontrast = 0;
        stimpara.onoffstep.greencontrast = 0.7;
        stimpara.onoffstep.bluecontrast = -0.7;
        stimpara.freqsweep.redcontrast = 0;
        stimpara.freqsweep.greencontrast = 0.7;
        stimpara.freqsweep.bluecontrast = -0.7;
        stimpara.contrastsweep.redcontrast = 0;
        stimpara.contrastsweep.greencontrast = 0.7;
        stimpara.contrastsweep.bluecontrast = -0.7;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        stimpara.coneisolating = true;
        
    case {'chirpnf30pf120frq300lf05hf8ctr300frq5cyanblack'}
        stimpara.color = true;
        stimpara.onoffstep.redcontrast = 0;
        stimpara.onoffstep.greencontrast = 0.7;
        stimpara.onoffstep.bluecontrast = 0.7;
        stimpara.freqsweep.redcontrast = 0;
        stimpara.freqsweep.greencontrast = 0.7;
        stimpara.freqsweep.bluecontrast = 0.7;
        stimpara.contrastsweep.redcontrast = 0;
        stimpara.contrastsweep.greencontrast = 0.7;
        stimpara.contrastsweep.bluecontrast = 0.7;
        [stimpara.redmeanintensity,stimpara.greenmeanintensity ,stimpara.bluemeanintensity] = deal(0,0.5,0.5);
        stimpara.coneisolating = true;
        
    case {'chirpnf180pf240frq480lf0hf8ctr480frq2rep15euler','opto_chirpnf180pf240frq480lf0hf8ctr480frq2rep15euler'}
        stimpara.Nrepeats = 15;
        stimpara.meanintensity = 0.5;
        stimpara.onoffstep.duration = 180;
        stimpara.onoffstep.preframes = 240;
        stimpara.onoffstep.contrast = 1;
        stimpara.freqsweep.preframes = 120;
        stimpara.freqsweep.duration = 480;
        stimpara.freqsweep.low_freq = 0.0;
        stimpara.freqsweep.high_freq = 8.0;
        stimpara.freqsweep.contrast = 1;
        stimpara.contrastsweep.preframes = 120;
        stimpara.contrastsweep.duration = 480;
        stimpara.contrastsweep.freq = 2.0;
        stimpara.contrastsweep.contrast = 1;
        stimpara.color = false;
        stimpara.coneisolating = false;
        stimpara.AfterEulerpaper = true;
        
    case 'chirpstimulus'
        stimpara.Nrepeats = 10;
        stimpara.onoffstep.preframes = 30;
        stimpara.onoffstep.contrast = 0.4;
        stimpara.freqsweep.preframes = 60;
        stimpara.freqsweep.contrast = 0.4;
        stimpara.contrastsweep.preframes = 60;
        stimpara.contrastsweep.freq = 10.0;
        stimpara.contrastsweep.contrast = 0.4;
        
    case 'chirpstimulus_lap4'
        stimpara.Nrepeats = 10;
        stimpara.onoffstep.preframes = 30;
        stimpara.onoffstep.contrast = 0.4;
        stimpara.freqsweep.preframes = 60;
        stimpara.freqsweep.contrast = 0.4;
        stimpara.contrastsweep.preframes = 60;
        stimpara.contrastsweep.freq = 10.0;
        stimpara.contrastsweep.contrast = 0.4;
        stimpara.pharmacology = true;
        stimpara.drugName = 'LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'chirpstimulus_cnqx'
        stimpara.Nrepeats = 10;
        stimpara.onoffstep.preframes = 30;
        stimpara.onoffstep.contrast = 0.4;
        stimpara.freqsweep.preframes = 60;
        stimpara.freqsweep.contrast = 0.4;
        stimpara.contrastsweep.preframes = 60;
        stimpara.contrastsweep.freq = 10.0;
        stimpara.contrastsweep.contrast = 0.4;
        stimpara.pharmacology = true;
        stimpara.drugName = 'CNQX and LAP4';
        stimpara.drugDose = '50 microMolar';
        
    case 'chirpstimulus_washout'
        stimpara.Nrepeats = 10;
        stimpara.onoffstep.preframes = 30;
        stimpara.onoffstep.contrast = 0.4;
        stimpara.freqsweep.preframes = 60;
        stimpara.freqsweep.contrast = 0.4;
        stimpara.contrastsweep.preframes = 60;
        stimpara.contrastsweep.freq = 10.0;
        stimpara.contrastsweep.contrast = 0.4;
        stimpara.pharmacology = true;
        stimpara.drugName = 'washout';
        stimpara.drugDose = 'Ames solution';
        
    case {'chirpstimulus_marmoset_15', 'chirpstimulus_marmoset_15r', 'chirpstim_marmoset_15r'...
            'chirpstimulus_marmoset_15r_bright79','chirpstimulus_marmoset_15r_bright91'...
            'chirpsimulus_marmoset_15r','chirp_marmoset_15r','chirpstimulus_marmoset_15r_after',...
            'chirpstimulus_marmoset_15r_win55212','chirpstimulus_marmoset_15r_before1h',...
            'chirpstimulus_marmoset_15repeat','chirpstimulus_15rp'}
        stimpara.Nrepeats = 15;
        stimpara.meanintensity = 0.5;
        stimpara.onoffstep.duration = 60;
        stimpara.onoffstep.preframes = 120;
        stimpara.onoffstep.contrast = 1;
        stimpara.freqsweep.preframes = 60;
        stimpara.freqsweep.duration = 480;
        stimpara.freqsweep.low_freq = 0.0;
        stimpara.freqsweep.high_freq = 15.0;
        stimpara.freqsweep.contrast = 1;
        stimpara.contrastsweep.preframes = 60;
        stimpara.contrastsweep.duration = 480;
        stimpara.contrastsweep.freq = 4.0;
        stimpara.contrastsweep.contrast = 1;
        stimpara.color = false;
        stimpara.coneisolating = false;
        stimpara.AfterEulerpaper = false;
        
        if strcmpi(desc,'chirpstimulus_marmoset_15r_bright79')
            stimpara.OLEDbrightness = 79;
        elseif strcmpi(desc,'chirpstimulus_marmoset_15r_bright91')
            stimpara.OLEDbrightness = 91;
        end
        
    case {'chirpstimulus_bkg0contrast1_75hz'
            'chirpstimulus_bkg0contrast1_75hz_40'
            'chirpstimulus_bkg0contrast1_75hz_0'}        
        stimpara.Nrepeats = 100;
        stimpara.meanintensity = 0.5;        
        stimpara.onoffstep.duration = 225;
        stimpara.onoffstep.preframes = 300;
        stimpara.onoffstep.contrast = 1;
        stimpara.onoffstep.preframesintensity = 0;        
        stimpara.freqsweep.preframes = 150;
        stimpara.freqsweep.duration = 600;
        stimpara.freqsweep.low_freq = 0.5;
        stimpara.freqsweep.high_freq = 8.0;
        stimpara.freqsweep.contrast = 1;        
        stimpara.contrastsweep.preframes = 225;
        stimpara.contrastsweep.duration = 600;
        stimpara.contrastsweep.freq = 2.0;
        stimpara.contrastsweep.contrast = 1;
        stimpara.color = false;
        stimpara.coneisolating = false;
        stimpara.AfterEulerpaper = true;  
end

newStimpara = stimpara;

end
