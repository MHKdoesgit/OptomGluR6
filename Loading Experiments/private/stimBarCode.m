

function [ newStimpara ] = stimBarCode (desc, stimpara)
%
%%% stimBarCode %%%
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
stimpara.stimulus = lower('barcode');
stimpara.preframes = 60;
stimpara.speed = 2;
stimpara.contrast = 1;
stimpara.maxperiod  = 1280;
stimpara.minperiod = 16;
stimpara.barcodelength = 1280;
stimpara.usevertical = true;
stimpara.meanintensity = 0.5;
stimpara.seed = -10000;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)    
    
    case 'barcodestimulus_plusvertical_120preframe'
        stimpara.preframes  = 1;
        stimpara.seed = -10074;
        
    case 'barcodestimulus_seed20000_for75hz'
        stimpara.preframes = 150;
        stimpara.seed = -20000;
        stimpara.contrast = 1;
        stimpara.usevertical = false;
end

newStimpara = stimpara;

end
