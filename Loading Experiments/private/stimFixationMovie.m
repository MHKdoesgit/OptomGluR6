

function [ newStimpara ] = stimFixationMovie (desc, stimpara)
%
%%% stimFixationMovie %%%
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


switch lower(desc)
    
    case {'fixationmovie_agergb_part1'}
        stimpara.stimulus = lower('fixationmovieAgeRGB_part1');
        stimpara.fixationpath = 'C:\NaturalStimulusFiles\age_rgb_part1\';
        stimpara.Nx = 1280;
        stimpara.Ny = 960;
        stimpara.frozenfixations = 126;
        stimpara.runningfixations = 800;
        
    case {'fixationmovie_agergb_part2'}
        stimpara.stimulus = lower('fixationmovieAgeRGB_part2');
        stimpara.fixationpath = 'C:\NaturalStimulusFiles\age_rgb_part2\';
        stimpara.Nx = 1280;
        stimpara.Ny = 960;
        stimpara.frozenfixations = 126;
        stimpara.runningfixations = 800;
        
    case {'fixationmovie_agegray_part1'}
        stimpara.stimulus = lower('fixationmovieAgeGRAY_part1');
        stimpara.fixationpath = 'C:\NaturalStimulusFiles\age_gray_part1\';
        stimpara.Nx = 1280;
        stimpara.Ny = 960;
        stimpara.frozenfixations = 112;
        stimpara.runningfixations = 800;
        
    case {'fixationmovie_agegray_part2'}
        stimpara.stimulus = lower('fixationmovieAgeGRAY_part2');
        stimpara.fixationpath = 'C:\NaturalStimulusFiles\age_gray_part2\';
        stimpara.Nx = 1280;
        stimpara.Ny = 960;
        stimpara.frozenfixations = 112;
        stimpara.runningfixations = 800;
        
    case {'fixationmovie_mouse_for75hz375run1875freeze',...
            'fixationmovie_mouse_for75hz375run1875freeze_40br',...
            'fixationmovie_mouse_for75hz375run1875freeze_0br'}
        
        stimpara.stimulus = lower('mouse_test_75Hz');
        stimpara.fixationpath = 'C:\NaturalStimulusFiles\mouse_test_75Hz\';
        stimpara.Nx = 1532;
        stimpara.Ny = 1024;
        stimpara.frozenfixations = 1875;
        stimpara.runningfixations = 375;
        stimpara.color = false;
   
end

newStimpara = stimpara;

end