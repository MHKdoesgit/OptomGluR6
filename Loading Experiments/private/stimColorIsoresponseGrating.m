
function [ newStimpara ] = stimColorIsoresponseGrating (desc, stimpara)
%
%%% stimColorIsoresponseGrating %%%
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
    
    case {'chromisorg24ang9cont1rep2to72cont60gr600std60per'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        
    case {'chromisorg24ang7cont1rep2to42cont60gr600std60per'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        stimpara.numrepeats  = 1;
        stimpara.numangles = 24;
        stimpara.numcontrasts = 7; 
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.42; % to 42% contrast      
end

newStimpara = stimpara;

end