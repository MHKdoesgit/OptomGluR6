

function [ newStimpara ] = stimPinkStripes (desc, stimpara)
%
%%% stimPinkStripes %%%
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
    
    case {'pinknoisestripes_0.1white'
            'pinkstripes0.1white'}
        stimpara.fracwhite = 0.1;
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 600;
        
    case {'pinknoisestripes_0.5white'
            'pinkstripes0.5white'}
        stimpara.fracwhite = 0.5;
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 600;
        
    case {'pinknoise10x10_2blinks_0.1white'
            'pinknoise10x10_2bl'
            'pinknoise10x10_2blinks'
            'pinknoise10x10_0.1white'
            'pinknoise0.1white'},
        stimpara.stixelwidth = 10;
        stimpara.stixelheight = 10;
        stimpara.fracwhite = 0.1;
        
    case 'temporalpinknoise',
        stimpara.stimulus = 'temporalcolorednoise';
        stimpara.temporalbeta = 1;
        stimpara.nodes = 300;
        stimpara.fracwhite = 0.1;
        stimpara.contrastwhite = 0.3;
        stimpara.contrastpink = 0.3;
        stimpara.meanintensity = 300;
end;

newStimpara = stimpara;
end

