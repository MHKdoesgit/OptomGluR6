

function [ newStimpara ] = stimOmittedStimulusResponse (desc, stimpara)
%
%%% stimOmittedStimulusResponse %%%
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
stimpara.stimulus = lower('omittedstimulusresponse');
stimpara.pulseframes = 2;
stimpara.pulseperiod = 5;
stimpara.preframes = 60;
stimpara.bkgintensity  = 0.5;
stimpara.pulseintensity = 0;
stimpara.seed = -1000;
stimpara.pulsenumbers = [4,8,16,32];
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)
    
    
    case 'omittedstimulusresponse_black_on_white'
        stimpara.bkgintensity  = 1;
        stimpara.pulseintensity = 0;
        
    case 'omittedstimulusresponse_white_on_black'
        stimpara.bkgintensity  = 0;
        stimpara.pulseintensity = 1;
end

newStimpara = stimpara;

end
