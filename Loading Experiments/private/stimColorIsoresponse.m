

function [ newStimpara ] = stimColorIsoresponse (desc, stimpara)
%
%%% stimColorIsoresponse %%%
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
    
    case {'chromisoresp36angle18cont3rep2to72cont30stim90pfr'}
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        
    case {'chromisoresp36angle10cont4rep2to42cont30stim90pfr'}
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.numrepeats  = 4;
        stimpara.numangles = 36;
        stimpara.numcontrasts = 10;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.42; % to 42% contrast
        
    case {'chromisoresp36angle15cont4rep2to72cont30stim90pfr'}
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.numrepeats  = 4;
        stimpara.numangles = 36;
        stimpara.numcontrasts = 15;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.72; % to 72% contrast
        
    case {'chromisoresp36angle15cont4rep2to72cont15stim45pfr'}
        stimpara.stimduration = 15;
        stimpara.preframes = 45;
        stimpara.numrepeats  = 4;
        stimpara.numangles = 36;
        stimpara.numcontrasts = 15;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.72; % to 72% contrast
        
    case {'chromisoresp36angle10cont4rep2to72cont15stim45pfr'}
        stimpara.stimduration = 15;
        stimpara.preframes = 45;
        stimpara.numrepeats  = 4;
        stimpara.numangles = 36;
        stimpara.numcontrasts = 10;
        stimpara.mincontrast = 0.02; % from 2% contrast
        stimpara.maxcontrast = 0.72; % to 72% contrast
end

newStimpara = stimpara;

end