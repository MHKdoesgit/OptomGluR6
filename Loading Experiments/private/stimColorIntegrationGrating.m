

function [ newStimpara ] = stimColorIntegrationGrating (desc, stimpara)
%
%%% stimColorIntegrationGrating %%%
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
    
    case {'colorintegrationgrating60gw600stimdur60period'
            'colorintegrationgratinggabazine60gw600stimdur60period'
            'colorintegrationgratingtpmpa60gw600stimdur60period'
            'colorintegrationgratinglap460gw600stimdur60period'
            'colorintegrationgratingstrychnine60gw600stimdur60period'
            'colorintegrationgratingwashout60gw600stimdur60period'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        
    case {'colorintegrationgrating30gw600stimdur60period'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 30;
        
    case {'colorintegrationgrating6gw600stimdur60period'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 6;
        
    case {'colorintegrationgrating2gw600stimdur60period'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 2;
        
    case {'colorintegrationgrating60gw600-40-4-40stimdur60period'}
        stimpara.stimduration = 600;
        stimpara.period = 60;
        stimpara.gratingwidth = 60;
        stimpara.contrastdiff = 0.04;
        stimpara.mincontrast = -0.4; % -40% contrast
        stimpara.maxcontrast = 0.4; % 40% contrast        
end

newStimpara = stimpara;

end

