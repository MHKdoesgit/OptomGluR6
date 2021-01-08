

function [ newStimpara ] = stimDriftingObjects (desc, stimpara)
%
%%% stimDriftingObjects %%%
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
    
    case 'driftingobjects'
        stimpara.seed = -10000;
        stimpara.stixelwidth = 1;
        stimpara.stixelheight = 1;
        stimpara.squarewidth = 800;
        stimpara.squareheight = 800;
        stimpara.radius = 14;
        stimpara.overlap = 0.2;
        stimpara.area = 303;
        stimpara.nobjects = ceil(stimpara.area/(2*stimpara.radius*(1-stimpara.overlap)));
        stimpara.distance = stimpara.squarewidth / stimpara.nobjects;
        stimpara.ndir = 8;
        stimpara.repeats = 5;
        stimpara.cycles = 2;
        stimpara.speed = 1.2;
        stimpara.meanintensity = 0.5;
        stimpara.contrast = 1;
        stimpara.preFrames = 100;
        stimpara.bkgcontrast = 0.3;
        stimpara.sequence = false;
        stimpara.eyemovements = false;
        stimpara.relativemotion = false;
        stimpara.bkgnoise = 0;
        stimpara.bkgstixel = 5;
        stimpara.bkgstep = 1;
        if stimpara.bkgnoise,   stimpara.bkgstixel = 10; end;
        if stimpara.relativemotion,
            stimpara.bkgnoise = 3;
            stimpara.bkgstixel = 40;
            stimpara.bkgstep = 2 * stimpara.speed;
        end;
end;

newStimpara = stimpara;

end

