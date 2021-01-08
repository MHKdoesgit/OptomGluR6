

function [ newStimpara ] = stimLinearRamp (desc, stimpara)
%
%%% stimLinearRamp %%%
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
    
    case 'linearramps_200period'
        stimpara.period = 200;
        
    case {'linearramps_period500', 'linearramps500period','linearramps_period500frames'}
        stimpara.period = 500;
        
    case {'linearramp_period600','linearramps600periodafterevenlonger'}
        stimpara.period = 600;
        
    case 'linearramps_period900'
        stimpara.period = 900;
end

newStimpara = stimpara;

end