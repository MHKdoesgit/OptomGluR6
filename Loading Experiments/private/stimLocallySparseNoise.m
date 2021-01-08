

function [ newStimpara ] = stimLocallySparseNoise (desc, stimpara)
%
%%% stimLocallySparseNoise %%%
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
stimpara.stimulus = lower('LocallySparseNoise');
stimpara.stimduration = 30;
stimpara.stixelwidth = 20;
stimpara.stixelheight = 20;
stimpara.gapwidth = 2;
stimpara.gapheight = 2;
stimpara.seedlocation = -1000;
stimpara.seedonoroffpixels = -10000;
stimpara.meanintensity = 0.5;
stimpara.contrast = 1;
stimpara.lmargin = 0;
stimpara.rmargin = 0;
stimpara.bmargin = 0;
stimpara.tmargin = 0;

switch lower(desc)
    
    case {'locallysparsenoise20x20gap2x2dur30circle'}
        stimpara.drawcircle = true;
        stimpara.drawannulus = false;
        stimpara.annulusstixel = 0;
        
    case {'locallysparsenoise20x20gap2x2dur60circle','opto_locallysparsenoise20x20gap2x2dur60circle'}
        stimpara.drawcircle = true;
        stimpara.drawannulus = false;
        stimpara.annulusstixel = 0;
        stimpara.stimduration = 60;
        
    case {'locallysparsenoise20x20gap2dur30','locallysparsenoise20x20gap2dur30frames',...
            'locallysparsenoise20x20gap2dur30f','locallysparsenoise_20x20gap2dur30f',...
            'locallysparsenoise20x20gap2dur30_am251_','locallysparsenoise_20x20gap2dur30_after',...
            'locallysparsenoise20x20gap2dur30f_control','locallysparsenoise20x20_gap2_dur30_win55212',...
            'locallysparsenoise20x20gap2dur30_after','locallysparsenoise_20x20gap2dur30',...
            'locallysparsenoise20x20_gap2dur30'}
        stimpara.drawcircle = false;
        stimpara.drawannulus = false;
        stimpara.annulusstixel = 0;
        stimpara.stimduration = 30;
        
    case {'locallysparsenoise20x20gap2x2surround40dur2annulus'}
        stimpara.drawcircle = false;
        stimpara.drawannulus = true;
        stimpara.annulusstixel = 40;
        stimpara.stimduration = 2;
        
    case {'locallysparsenoise20x20gap1x1surround40dur60annulus'}
        stimpara.drawcircle = false;
        stimpara.drawannulus = true;
        stimpara.annulusstixel = 40;
        stimpara.stimduration = 60;
        stimpara.gapwidth = 1;
        stimpara.gapheight = 1;
        
    case {'coloropponentcentersurround20cent60surr60stimdur'}
        stimpara.stimulus = 'coloropponentcentersurround';
        stimpara.stimduration = 60;
        stimpara.centerdimeter = 20;
        stimpara.surrounddiameter = 60;
        stimpara.gapwidth = 1;
        stimpara.gapheight = 1;
        stimpara.contrast = 0.75;
        stimpara.redmean = 0.0;
        stimpara.greenmean = 0.5;
        stimpara.bluemean = 0.5;
        stimpara.seedstimseq = -2000;
        stimpara.coneisolation = true;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight','meanintensity'});
        stimpara.sparsegridsize = stimpara.surrounddiameter;
        
    case {'coloropponentcentersurround20cent60surr60stimdur3x3gap'}
        stimpara.stimulus = 'coloropponentcentersurround';
        stimpara.stimduration = 60;
        stimpara.centerdimeter = 20;
        stimpara.surrounddiameter = 60;
        stimpara.gapwidth = 3;
        stimpara.gapheight = 3;
        stimpara.contrast = 0.75;
        stimpara.redmean = 0.0;
        stimpara.greenmean = 0.5;
        stimpara.bluemean = 0.5;
        stimpara.seedstimseq = -2000;
        stimpara.coneisolation = true;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight','meanintensity'});
        stimpara.sparsegridsize = stimpara.centerdimeter;
        
    case 'coloropponentcentersurround25cent75surr30stimdur4x4gap'
        stimpara.stimulus = 'coloropponentcentersurround';
        stimpara.stimduration = 30;
        stimpara.centerdimeter = 25;
        stimpara.surrounddiameter = 75;
        stimpara.gapwidth = 4;
        stimpara.gapheight = 4;
        stimpara.contrast = 0.75;
        stimpara.redmean = 0.0;
        stimpara.greenmean = 0.5;
        stimpara.bluemean = 0.5;
        stimpara.seedstimseq = -2000;
        stimpara.coneisolation = true;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight','meanintensity'});
        stimpara.sparsegridsize = stimpara.centerdimeter;
        
    case 'coloropponentcentersurround20cent60surr15stimdur'
        stimpara.stimulus = 'coloropponentcentersurround';
        stimpara.stimduration = 15;
        stimpara.centerdimeter = 20;
        stimpara.surrounddiameter = 60;
        stimpara.gapwidth = 3;
        stimpara.gapheight = 3;
        stimpara.contrast = 0.75;
        stimpara.redmean = 0.0;
        stimpara.greenmean = 0.5;
        stimpara.bluemean = 0.5;
        stimpara.seedstimseq = -2000;
        stimpara.coneisolation = true;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight','meanintensity'});
        stimpara.sparsegridsize = stimpara.centerdimeter;
        
    case {'sparsebar50w6h4angle1gap60dur'}
        stimpara.stimulus = 'sparsebar';
        stimpara.barwidth = 50;
        stimpara.barheight = 6;
        stimpara.numangles = 8;
        stimpara.gapwidth = 1;
        stimpara.gapheight = 1;
        stimpara.seedstimseq = -2000;
        stimpara.stimduration = 60;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight'});
        
    case {'sparsespot8162550rad1gap60dur'}
        stimpara.stimulus = 'sparsespot';
        stimpara.spotdiameters = [8 16 25 50];
        stimpara.gapwidth = 1;
        stimpara.gapheight = 1;
        stimpara.seedstimseq = -2000;
        stimpara.stimduration = 60;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight'});
        stimpara.sparsegridsize = max(stimpara.spotdiameters);
        
    case {'sparsespot8162550rad7gap60dur'}
        stimpara.stimulus = 'sparsespot';
        stimpara.spotdiameters = [8 16 25 50];
        stimpara.gapwidth = 7;
        stimpara.gapheight = 7;
        stimpara.seedstimseq = -2000;
        stimpara.stimduration = 60;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight'});
        stimpara.sparsegridsize = min(stimpara.spotdiameters);
        
    case {'sparsespot162550rad1gap60dur','sparsespot162550rad7gap60dur','sparsespot162550rad4gap60dur'}
        stimpara.stimulus = 'sparsespot';
        stimpara.spotdiameters = [16 25 50];
        stimpara.gapwidth = 4;
        stimpara.gapheight = 4;
        stimpara.seedstimseq = -2000;
        stimpara.stimduration = 60;
        stimpara.changecontperframe = false;
        stimpara = rmfield(stimpara,{'stixelwidth', 'stixelheight'});
        stimpara.sparsegridsize = min(stimpara.spotdiameters);
        
end

newStimpara = stimpara;

end
