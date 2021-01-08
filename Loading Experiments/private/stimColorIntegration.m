

function [ newStimpara ] = stimColorIntegration (desc, stimpara)
%
%%% stimColorIntegration %%%
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
    
    case {'colorintegration-20step2to20cont30stim120prefr'
            'colorintegration-20step2to20cont30stimdur120prefr'}
        stimpara.stimduration = 30;
        stimpara.preframes = 120;
        
    case {'colorintegration-20step2to20cont30stimdur90prefr'
            'colorintegrationbeforedrug-20-2-20cont30stim90prefr'
            'colorintegrationgabazine-20-2-20cont30stim90prefr'
            'colorintegrationtpmpa-20-2-20cont30stim90prefr'
            'colorintegrationlap4-20-2-20cont30stim90prefr'
            'colorintegrationstrychnine-20-2-20cont30stim90prefr'
            'colorintegrationwashout-20-2-20cont30stim90prefr'
            'colorintegration-20step2to20cont30stim90pfr_highphotopic'            
            'colorintegration-20step2to20cont30stimdur90prefr_tpmpa'
            'colorintegrationhepes-20-2-20cont30stim90prefr'
            'colorintegrationhepes-20step2to20cont30stimdur90prefr'}
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        
    case 'colorintegration-40step4to40cont30stimdur90prefr'
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.contrastdiff = 0.04;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4;  % 40% contrast
        
    case 'colorintegration-60step6to60cont30stimdur90prefr'
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.contrastdiff = 0.06;
        stimpara.mincontrast = -0.6; % -60% contrast
        stimpara.maxcontrast = 0.6;  % 60% contrast
        
    case 'colorintegration-60step2to60cont30stimdur90prefr'
        stimpara.stimduration = 30;
        stimpara.preframes = 90;
        stimpara.contrastdiff = 0.02; % step of 2 for more combinations
        stimpara.mincontrast = -0.6; % -60% contrast
        stimpara.maxcontrast = 0.6;  % 60% contrast
        
    case {'sparsecolorintegration20x20gap2x2dur30circle'}
        stimpara.stimulus = lower('SparseColorIntegration');
        stimpara.stimduration = 30;
        stimpara = rmfield(stimpara, {'preframes','seed'});
        stimpara.stixelwidth = 20;
        stimpara.stixelheight = 20;
        stimpara.gapwidth = 2;
        stimpara.gapheight = 2;
        stimpara.seedlocation = -1000;
        stimpara.seedonoroffpixels = -10000;
        stimpara.seedcolorcontrast = -2000;
        stimpara.drawcircle = true;
        stimpara.drawannulus = false;
        
    case 'sparsecolorintegration20x20gap2x2dur15circle'
        stimpara.stimulus = lower('SparseColorIntegration');
        stimpara.stimduration = 15;
        stimpara = rmfield(stimpara, {'preframes','seed'});
        stimpara.stixelwidth = 20;
        stimpara.stixelheight = 20;
        stimpara.gapwidth = 2;
        stimpara.gapheight = 2;
        stimpara.seedlocation = -1000;
        stimpara.seedonoroffpixels = -10000;
        stimpara.seedcolorcontrast = -2000;
        stimpara.drawcircle = true;
        stimpara.drawannulus = false;
        
    case 'sparsecolorintegration-40-4-40cont20x20gap2x2dur15circle'
        stimpara.stimulus = lower('SparseColorIntegration');
        stimpara.stimduration = 15;
        stimpara = rmfield(stimpara, {'preframes','seed'});
        stimpara.contrastdiff = 0.04;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4; % 40% contrast
        stimpara.stixelwidth = 20;
        stimpara.stixelheight = 20;
        stimpara.gapwidth = 2;
        stimpara.gapheight = 2;
        stimpara.seedlocation = -1000;
        stimpara.seedonoroffpixels = -10000;
        stimpara.seedcolorcontrast = -2000;
        stimpara.drawcircle = true;
        stimpara.drawannulus = false;    
end

newStimpara = stimpara;

end

