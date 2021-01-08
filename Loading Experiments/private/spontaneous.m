

function [ newStimpara ] = spontaneous (desc, stimpara)
%
%%% spontaneous %%%
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
    
    case {'spontaneous','spontaneous_gray','spontaneous_activity_gray','spontaneous_activity_grey'...
            'spontanous_activity','after2_spontact_gray','spontactlight_adapt'...
            'spontaneous_activity_grey_beforedrug','spontactivity_gray_afterlongtime'...
            'spontaneous_activity_grey_washin','spontactivity_gray_washin'...
            'spontaneous_activity_grey_washout','spontaneousactivity_gray_washout' ...
            'opto_spontaneous_activity_grey_bleaching'...
            'opto_spontaneous_activity_grey_ndf3','spontactivity_gray_adaptation'...
            'spontact_gray', 'spontactivity_gray_long','goingbacktogray','turning_oled_on'...
            'spontactivity_gray''spontaneousactivity_gray','steppingtogray'...
            'spontactivitygray_long','spontactivitygray','goingtogray',...
            'spontactivity_gray_washout','spontactivity_gray_after','spontaneousactivity_darkadapt',...
            'spontaneousactivity_darkadapt2'}
        
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        
    case {'opto_spontaneous_activity_white_bleaching', 'spontaneous_light','spontactivity_light',...
            'spontact_light'}
        stimpara.stimulus = 'spontaneous_white';
        stimpara.contrast = 1;
        
    case {'spontaneous_black','spontactivity_dark', 'spontaneousactivity_dark','spontactivity_dark1'...
            'spontaneous_activity_black','spontactivity_darkadapt','spontaneous_dark','turningmonitoron'...
            'spontactivity_darkadapt1','spontactivity_darkadapt2','spontaneousdark','dark'...
            'spontactdark_duringadapt'}
        
        stimpara.stimulus = 'spontaneous_black';
        stimpara.contrast = 0;
        
    case 'spontactivity_long_bright43'
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        stimpara.OLEDbrightness = 43;
        
    case 'spontactivity_bright55'
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        stimpara.OLEDbrightness = 55;
        
    case    'spontactivity_bright67'
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        stimpara.OLEDbrightness = 67;
        
    case    'spontactivity_bright79'
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        stimpara.OLEDbrightness = 79;
        
    case    'spontactivity_bright91'
        stimpara.stimulus = 'spontaneous_gray';
        stimpara.contrast = 0.5;
        stimpara.OLEDbrightness = 91;
end

newStimpara = stimpara;

end