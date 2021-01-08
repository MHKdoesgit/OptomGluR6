

function [ newStimpara ] = stimConeIsolationTest (desc, stimpara)
%
%%% stimConeIsolationTest %%%
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
    
    case 'coneisotest_1_0gmg0.5_30frames120preframes'
        stimpara.greencontrast = 0;
        stimpara.greenmeanintensity = 0.5;
        
    case 'coneisotest_2_10gmg0.5_30frames120preframes'
        stimpara.greencontrast = -0.1;
        stimpara.greenmeanintensity = 0.5;
        
    case 'coneisotest_3_20gmg0.5_30frames120preframes'
        stimpara.greencontrast = -0.2;
        stimpara.greenmeanintensity = 0.5;
        
    case 'coneisotest_4_30gmg0.5_30frames120preframes'
        stimpara.greencontrast = -0.3;
        stimpara.greenmeanintensity = 0.5;
        
    case 'coneisotest_5_40gmg0.5_30frames120preframes'
        stimpara.greencontrast = -0.4;
        stimpara.greenmeanintensity = 0.5;
        
    case 'coneisotest_6_0gmg0.6_30frames120preframes'
        stimpara.greencontrast = 0;
        stimpara.greenmeanintensity = 0.6;
        
    case 'coneisotest_7_10gmg0.6_30frames120preframes'
        stimpara.greencontrast = -0.1;
        stimpara.greenmeanintensity = 0.6;
        
    case 'coneisotest_8_20gmg0.6_30frames120preframes'
        stimpara.greencontrast = -0.2;
        stimpara.greenmeanintensity = 0.6;
        
    case 'coneisotest_9_30gmg0.6_30frames120preframes'
        stimpara.greencontrast = -0.3;
        stimpara.greenmeanintensity = 0.6;
        
    case 'coneisotest_10_40gmg0.6_30frames120preframes'
        stimpara.greencontrast = -0.4;
        stimpara.greenmeanintensity = 0.6;
        
    case 'coneisotest_10to40grcont_06newmean_70uvcont'
        stimpara.stimduration = 30;
        stimpara.preframeduration = 90;
        stimpara.redcontrast = 0;
        stimpara.bluecontrast = 0.7;
        stimpara.greencontrastdiff = 0.1;
        stimpara.greenmincontrast = 0.0;
        stimpara.greenmaxcontrast = 0.4;
        stimpara.newgreenmean = 0.6;
        stimpara.seed = -1000;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
        
end

newStimpara = stimpara;

end
