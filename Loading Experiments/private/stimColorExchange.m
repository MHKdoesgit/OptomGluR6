

function [ newStimpara ] = stimColorExchange (desc, stimpara)
%
%%% stimColorExchange %%%
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
    
    case {'colorexchange-40step20to40with20fix30stimdur90pfr'}
        stimpara.contrastdiff = 0.20;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4; % 40% contrast
        stimpara.fixedcontrast = 0.20; % 10% fixed contrast
        
    case {'colorexchange-40step10to40with10fix30stimdur90pfr'}
        stimpara.contrastdiff = 0.10;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4; % 40% contrast
        stimpara.fixedcontrast = 0.10; % 10% fixed contrast      
end

newStimpara = stimpara;

end

