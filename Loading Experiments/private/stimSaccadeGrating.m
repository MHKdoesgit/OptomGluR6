

function [ newStimpara ] = stimSaccadeGrating (desc, stimpara)
%
%%% stimSaccadeGrating %%%
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
    
    case {'saccadegrating_18pixshift2_100ms_for60hz'
            'saccadegratings_18pixshift2_100ms'}
        stimpara.fixationframes = 48;
        stimpara.saccadeframes = 6;
        stimpara.barwidth = 18;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
        
    case 'saccadegrating_18pix_400fixshift2_100ms_for60hz'
        stimpara.fixationframes = 400;
        stimpara.saccadeframes = 6;
        stimpara.barwidth = 18;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
        
    case {'saccadegrating_16pixshift2_100ms','opto_saccadegrating_16pixshift2_100ms'}
        stimpara.fixationframes = 48;
        stimpara.saccadeframes = 6;
        stimpara.barwidth = 16;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
        
    case 'saccadegrating_12pixshift2_67ms'
        stimpara.fixationframes = 32;
        stimpara.saccadeframes = 4;
        stimpara.barwidth = 12;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
        
    case {'saccadegrating12pixshift2_50ms', 'saccadegrating_12pixshift2_50',...
            'saccadegrating_12pixshift2_50ms','saccadegratings_12pixshift2_50ms'}
        stimpara.fixationframes = 21;
        stimpara.saccadeframes = 3;
        stimpara.barwidth = 12;
        stimpara.cosineprofile = false;
        stimpara.averageshift = 2;
end

newStimpara = stimpara;

end

