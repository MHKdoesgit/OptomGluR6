

function [ newStimpara ] = stimObjectsMovingBackground (desc, stimpara)
%
%%% ObjectsMovingBackground %%%
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
stimpara.stimulus = lower('objectsmovingbackground');
stimpara.nblinks = 2;
stimpara.seed = -1000;
stimpara.cycles = 1;
stimpara.preframes = 200;
stimpara.stimFrames = 108000;
stimpara.addXoffset = 0;
stimpara.addYoffset = 0;
stimpara.stixelwidth = 1;
stimpara.stixelheight = 1;
stimpara.nobjects = -1;
stimpara.radius = 16;
stimpara.distance = 5 * 16;
stimpara.overlap = 0;
stimpara.objectdist = 1;
stimpara.objectmotion = 0;
stimpara.resetforce = 0;
stimpara.speed = 1;
stimpara.eyemovements = true;
stimpara.gausssteps = true;
stimpara.bgnoise = 4;
stimpara.bgcontrast = 0.3;
stimpara.stepsize = 2;
stimpara.stepsizemultiplier = 1;
stimpara.bgstixel = 5;
stimpara.filterstdv = 5;
stimpara.filtermultiplier = 1;
stimpara.pinknoiseslope = 1;
stimpara.smoothtrajectory = 0;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.imagePath = 'Y:\\Norma\\OpenGL_images\\Test_0001.raw';
stimpara.imageWidth = 1000;
stimpara.imageHeight = 900;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)
    
    case {'omb_bg1x1corr4stdvc150_gsteps3stdv','omb_bg1x1corr4stdv_c150_gsteps3stdv_15min'}
        stimpara.radius = 0;
        stimpara.stimFrames = 54000;
        stimpara.gausssteps = true;
        stimpara.stepsize = 3;
        stimpara.bgnoise = 4;
        stimpara.bgcontrast = 1.5;
        stimpara.bgstixel = 1;
        stimpara.filterstdv = 4;
        
    case {'omb_bg4x4corr4stdvc150_gsteps3stdv','omb_bg4x4corr4stdv_c150_gsteps3stdv_15min',...
            'omb_bg4x4corr8stdv_c150_gsteps3stdv'}
        stimpara.radius = 0;
        stimpara.stimFrames = 54000;
        stimpara.gausssteps = true;
        stimpara.stepsize = 3;
        stimpara.bgnoise = 4;
        stimpara.bgcontrast = 1.5;
        stimpara.bgstixel = 4;
        stimpara.filterstdv = 8;
end

newStimpara = stimpara;

end
