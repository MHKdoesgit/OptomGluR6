


function [ newStimpara ] = stimSpatialFrequencyContrastSensitivity (desc, stimpara)
%
%%% stimSpatialFrequencyContrastSensitivity %%%
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


switch lower(desc)
    
    case {'spatialfrequencysensitivity','spatialfrequencysensitivity_bright55'...
            'spatialfrequencysensitivity_bright67'...
            'spatialfrequencysensitivity_bright79'...
            'spatialfrequsensitivity_bright91'}
        stimpara.stimulus = lower('spatialfrequencysensitivity');
        stimpara.nrepeats = 1;
        stimpara.nframes = 600;
        stimpara.preframes = 120;
        stimpara.postframes = 60;
        stimpara.contrast = 1;
        stimpara.meanintensity = 0.5;
        stimpara.radius = 10;
        stimpara.temporalPeriod = 20;
        stimpara.centerX = 0;
        stimpara.centerY = 0;
        stimpara.untriggeredMode = true;
        stimpara.readFromFile = false;
        stimpara.circleMask  = false;
        stimpara.spatialperiod = [1730 173 109.5 68.9 43.5 27.4 17.3 10.95 6.89 4.35];
        
        if strcmpi(desc,'spatialfrequencysensitivity_bright55')
            stimpara.OLEDbrightness = 55;
        elseif strcmpi(desc,'spatialfrequencysensitivity_bright67')
            stimpara.OLEDbrightness = 67;
        elseif strcmpi(desc,'spatialfrequencysensitivity_bright79')
            stimpara.OLEDbrightness = 79;
        elseif strcmpi(desc,'spatialfrequsensitivity_bright91')
            stimpara.OLEDbrightness = 91;
        end
        
    case {'contrastsensitivitywithdriftinggratings'
            'contrastsensitivitywithdriftinggrating'
            'contrastsensitivity_with_drifting_gratings'
            'contrastsensitivity_withdriftinggratings'
            'contrastsensitivity_with_driftinggratings'
            'contrastsensitivtywithdriftinggratings'
            'contrastsensitivtywithdriftinggrating'
            'contrastsensitivitywithdriftingratings_control'
            'contrastsensitivitywithdriftingratings_am251'
            'contrastsensitivitywithdriftinggratings_control'
            'contrastsensitivitywithdriftinggratings_win55212'
            'contrastsensitivitywithdriftinggratings_after'
            'contrastsensitivitydriftinggrating'}
        stimpara.stimulus = lower('contrastsensitivitywithdriftinggratings');
        stimpara.preframes = 30;
        stimpara.nframes = 300;
        stimpara.postframes = 30;
        stimpara.nrepeats = 1;
        stimpara.meanintensity = 0.5;
        stimpara.temporalPeriod = 20;
        stimpara.spatialperiod = [1024 512 256 128 64 32 16 8 4];
        stimpara.contrastset = [1 2 4 8 16 32];
end

newStimpara = stimpara;

end
