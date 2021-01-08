

function [ newStimpara ] = stimContrastStepLadder (desc, stimpara)
%
%%% stimContrastStepLadder %%%
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
    
    case {'contraststeps'
            'contraststepladder'}
        stimpara.nframes = 12;
        stimpara.preframes = 108;
        stimpara.contrasts = [0, 0.01, 0.02,0.03, 0.05, 0.08, 0.12, 0.18, 0.24]; % From 20/02/2013 on
        disp('Using contraststepladder contrasts from 20/02/2013...');
        stimpara.colormode = 'black & white';
        
    case 'contraststepladder_12fr_180preframes_white'
        stimpara.color = false;
        stimpara.colormode = 'white';
        stimpara.coneisolating = false;
        
    case 'contraststepladder_12fr_180preframes_blue'
        stimpara.color = true;
        stimpara.colormode = 'uv';
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'contraststepladder_12fr_180preframes_green'
        stimpara.color = true;
        stimpara.colormode = 'green';
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'contraststepladder_12fr_180preframes_cyan'
        stimpara.color = true;
        stimpara.colormode = 'cyan';
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'contraststepladderhighercontrasts'
        stimpara.nframes = 12;
        stimpara.preframes = 108;
        stimpara.contrasts = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.16 0.23 0.32 0.45 0.64];
        stimpara.colormode = 'black & white';
        
    case 'contraststepladder_12fr_108preframes_cyan'
        stimpara.nframes = 12;
        stimpara.preframes = 108;
        stimpara.color = true;
        stimpara.colormode = 'cyan';
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'contraststepladder_12fr_108preframes_highcont'
        stimpara.nframes = 12;
        stimpara.preframes = 108;
        stimpara.contrasts = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.16 0.23 0.32 0.45 0.64];
        stimpara.color = true;
        stimpara.colormode = 'cyan';
        stimpara.coneisolating = true;
        stimpara.usered = false;
        stimpara.usegreen = true;
        stimpara.useblue = true;
        stimpara.redmeanintensity = 0;
        stimpara.greenmeanintensity = 0.5;
        stimpara.bluemeanintensity = 0.5;
        
    case 'contraststepladder12preframes48levels8'
        stimpara.nframes = 12;
        stimpara.preframes = 48;
        stimpara.contrasts = [0.01 0.02 0.03 0.05 0.1 0.2 0.6 1.0];
        stimpara.colormode = 'black & white';
        
    case {'contraststepladder_12fr_108pfr_1234681624325064contrasts'
            'opto_contraststepladder_12fr_108pfr_1234681624325064c'}
        stimpara.nframes = 12;
        stimpara.preframes = 108;
        stimpara.contrasts = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.16 0.24 0.32 0.5 0.64];
        stimpara.colormode = 'black & white';
        
    case {'contraststepladder_12fr_288pfr_2481624325064contrasts'}
        stimpara.nframes = 12;
        stimpara.preframes = 288;
        stimpara.contrasts = [0 0.02 0.04 0.08 0.16 0.24 0.32 0.5 0.64];
        stimpara.colormode = 'black & white';
        
    case {'contraststepladder12preframe48levels8','contraststepladder12preframes48levels8n'...
            'cointraststepladder12preframes48levels8n','contraststepladder12preframes48levels8_control',...
            'contraststepsladder12preframes48levels8_am251','contraststepladder12preframes48levels8_after',...
            'contraststepladder12preframes48levels8_win55212','contraststepladdere12preframes48levels8_after',...
            'contraststepsladder12preframes48levels8','contraststepladder12preframes48_8levels',...
            'constraststepladder12preframes48levels8'}
        stimpara.nframes = 12;
        stimpara.preframes = 48;
        stimpara.contrasts = [0 0.01 0.02 0.03 0.05 0.1 0.2 0.6 1];
        stimpara.colormode = 'black & white';
        
    case 'contrastladder_30blinks_60preframes'
        stimpara.nframes = 30;
        stimpara.preframes = 60;
        stimpara.contrasts = [0 0.01 0.02 0.03 0.05 0.1 0.2 0.6 1];
        stimpara.colormode = 'black & white';
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
        
end

if strcmpi(stimpara.colormode,'black & white')
    stimpara = rmfield(stimpara,{'color','coneisolating','usered','usegreen','useblue',...
        'redmeanintensity','greenmeanintensity','bluemeanintensity'});
    stimpara.color = false;
end

newStimpara = stimpara;

end
