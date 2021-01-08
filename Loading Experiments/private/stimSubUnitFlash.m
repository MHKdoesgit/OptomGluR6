

function [ newStimpara ] = stimSubUnitFlash (desc, stimpara)
%
%%% stimSubUnitFlash %%%
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
stimpara.stimulus = lower('subunitflash');
stimpara.stixelwidth = 10;
stimpara.stixelheight = 10;
stimpara.emptyspace = 0;
stimpara.emptyspacew = 0;
stimpara.emptyspaceh = 0;
stimpara.nangles = 8;
stimpara.meanintensity = 0.5;
stimpara.ncontrasts = 10;
stimpara.mincontrast = 0;
stimpara.maxcontrast = 1;
stimpara.seed = -10000;
stimpara.stimduration = 120;
stimpara.flashstart = 30;
stimpara.flashstop = 42;
stimpara.nrepeats = 10;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)
    
    case {'subunitflash588_7lin'
            'subunitflash7_588_7lin'
            'subunitflash7_588lin24x10p0s2'
            'subunitflash21_196lin24x10p0s2'
            'subunitflash7_588log5p1s2'
            'subunitflash14_588lin20x12p0s2'
            'subunitflash7_588lin20x12p0s2'}
        
        stimpara.flashorder = 'michael';
        
    case {'suf_5pxspace_0.03_0.63_6cont_15x15'
            'suf_5pxspace_0.03_0.63_6cont_16angles_15x15'
            'suf_5pxspace_0.03_0.63_12cont_15x15'
            'suf_2pxspace_0.03_0.63_12con_16ang_15x15'
            'suf_nospace_0.03_0.63_12contrasts_16angles_15x15'
            'suf_5pxspace_0.03_0.63_6cont_7x7'
            'suf_2pxspace_0.03_0.63_12con_16ang_7x7'
            'suf_nospace_0.03_0.63_6cont_7x7'
            'suf_nospace_0.03_0.63_12contrasts_7x7'
            'suf_10px_0.03_0.63_12con_16ang_10x10'
            'suf_8px_0.03_0.63_12con_16ang_8x8'
            'suf_8px_0.03_0.63_12con_16ang_5x5'}
        
        stimpara.flashorder = 'Fernando';
        
        
    case {'subunitflash14x14ang24con10rep4dur200ms'
            'subunitflash14x14ang24con0rep4dur200ms'
            'subunitflash14x14ang24con10rep4dur200ms_500ms'}
        
        stimpara.stixelwidth = 14;
        stimpara.stixelheight = 14;
        stimpara.emptyspace = 0;
        stimpara.mincontrast = 0.03;
        stimpara.maxcontrast = 1;
        stimpara.nangles = 24;
        stimpara.stimduration = 60;
        stimpara.flashstart = 8;
        stimpara.flashstop = 19;
        stimpara.nrepeats = 4;
end

newStimpara = stimpara;

end
