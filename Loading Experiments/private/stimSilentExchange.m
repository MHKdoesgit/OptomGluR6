

function [ newStimpara ] = stimSilentExchange (desc, stimpara)
%
%%% stimSilentExchange %%%
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
    
    case 'silentexchange_1_20g0uv_30frames120preframes'        
        stimpara.greencontrast = 0.2;
        stimpara.bluecontrast = 0;
        
    case 'silentexchange_2_18g-2uv_30frames120preframes'        
        stimpara.greencontrast = 0.18;
        stimpara.bluecontrast = -0.02;
        
    case 'silentexchange_3_16g-4uv_30frames120preframes'        
        stimpara.greencontrast = 0.16;
        stimpara.bluecontrast = -0.04;
        
    case 'silentexchange_4_14g-6uv_30frames120preframes'        
        stimpara.greencontrast = 0.14;
        stimpara.bluecontrast = -0.06;
        
    case 'silentexchange_5_12g-8uv_30frames120preframes'        
        stimpara.greencontrast = 0.12;
        stimpara.bluecontrast = -0.08;
        
    case 'silentexchange_6_20uv0g_30frames120preframes'        
        stimpara.greencontrast = 0;
        stimpara.bluecontrast = 0.2;
        
    case 'silentexchange_7_18uv-2g_30frames120preframes'        
        stimpara.greencontrast = -0.02;
        stimpara.bluecontrast = 0.18;
        
    case 'silentexchange_8_16uv-4g_30frames120preframes'        
        stimpara.greencontrast = -0.04;
        stimpara.bluecontrast = 0.16;
        
    case 'silentexchange_9_14uv-6g_30frames120preframes'        
        stimpara.greencontrast = -0.06;
        stimpara.bluecontrast = 0.14;
        
    case 'silentexchange_10_12uv-8g_30frames120preframes'        
        stimpara.greencontrast = -0.08;
        stimpara.bluecontrast = 0.12;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
        
end

newStimpara = stimpara;

end


