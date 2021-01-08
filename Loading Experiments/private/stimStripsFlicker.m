

function [ newStimpara ] = stimStripsFlicker (desc, stimpara)
%
%%% stimStripsFlicker %%%
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
    
    case 'whitenoisestripes'
        stimpara.seed = -10000;
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 600;
        
    case 'whitenoisestripes_5x600'
        stimpara.seed = -10000;
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 600;
        
    case 'stripeflicker_4x600_bw_4blinks'
        stimpara.seed = -1000;
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 480;
        stimpara.nblinks = 4;
        
    case {'stripeflicker2x600bw1blink_control'
            'stripeflicker2x600bw1blink_am251_'
            'stripeflicker2x600bw1blink_after'
            'stripeflicker2x600bw1blink_win55212'
            'stripeflicker2x600bw1blink'}
        stimpara.seed = -1000;
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 600;
        stimpara.nblinks = 1;
        
    case {'stripeflicker5x600bw1blink'}
        stimpara.seed = -1000;
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 600;
        stimpara.nblinks = 1;
        
end

newStimpara = stimpara;

end
