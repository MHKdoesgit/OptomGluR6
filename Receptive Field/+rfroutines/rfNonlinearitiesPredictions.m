

function res = rfNonlinearitiesPredictions(res, runningbin, frozenbin, stimPara, varargin)

switch stimPara.colormode
    case 'monochromatic'
        res = nlpredmonochromatic(res, runningbin, frozenbin, stimPara);
        
    case 'dichromatic'
        res = nlpreddichromatic(res, runningbin, frozenbin, stimPara);
        
    case 'trichromatic'
        res = nlpredtrichromatic(res, runningbin, frozenbin, stimPara);
        
end

end

function res = nlpredmonochromatic(res, runningbin, frozenbin, stimPara)

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;
Ncells = size(clus,1);
Ntrials = size(res.trialRates,3);

fprintf('Calculating generator signals... '); tic;

temporalComponents = temporalComponents./sqrt(sum(temporalComponents.^2,2));
lrspacepredict  = reshape(flip(spatialComponents, 2), Ncells, Ny*Nx);
lrspacepredict  = lrspacepredict./sqrt(sum(lrspacepredict.^2,2));

% calculate low-rank generator signals
[lrgenerators, ~] = rf.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    lrspacepredict, temporalComponents, stimPara.seed);

modeltcomps  = modeltcomps./sqrt(sum(modeltcomps.^2,2));
modelspredict  = reshape(flip(modelscomps, 2), Ncells, Ny*Nx);
modelspredict  = modelspredict./sqrt(sum(modelspredict.^2,2));

% calculate model generator signals
[modelgenerators, spikes] = rf.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    modelspredict, modeltcomps, stimPara.seed);

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
fprintf('Extracting nonlinearities... '); tic;

nlncentslr    = NaN(Ncells, stimPara.nonlinBinN, 'single');
nlnvalslr     = NaN(Ncells, stimPara.nonlinBinN, 'single');
nlncentsmodel = NaN(Ncells, stimPara.nonlinBinN, 'single');
nlnvalsmodel  = NaN(Ncells, stimPara.nonlinBinN, 'single');

for icell = 1:Ncells
    
    cellgenslr = lrgenerators(icell,:)';
    if isnan(sum(cellgenslr(:))); continue; end
    [staVals, staCents,~] = rf.getNonlinearity(cellgenslr, spikes(icell,:),...
        stimPara.nonlinBinN,stimPara.Nblinks/stimPara.refreshrate);
    nlncentslr(icell,:) = staCents;
    nlnvalslr(icell, :) = staVals;
    
    cellgensmodel = modelgenerators(icell,:)';
    if isnan(sum(cellgensmodel)); continue; end
    [staVals,staCents,~] = rf.getNonlinearity(cellgensmodel, spikes(icell,:),...
        stimPara.nonlinBinN,stimPara.Nblinks/stimPara.refreshrate);
    nlncentsmodel(icell,:) = staCents;
    nlnvalsmodel(icell,:) = staVals;
end

res.nlncentslr    = nlncentslr;             res.nlnvalslr           = nlnvalslr;
res.nlncentsmodel = nlncentsmodel;          res.nlnvalsmodel        = nlnvalsmodel;

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
fprintf('Calculating predictions and performances... '); tic;

%Generate stimulus
[frozenstimulus, ~] = rf.ran1bool(stimPara.secondseed, stimPara.FrozenFrames*Nx*Ny);
frozenstimulus      = reshape(frozenstimulus, Ny*Nx, stimPara.FrozenFrames);
frozenstimulus      = 2 * single(frozenstimulus) - 1; % transform to contrast

% filter with spatial components
filteredstim_model = modelspredict  * frozenstimulus;
filteredstim_lr    = lrspacepredict * frozenstimulus;

trialrates = frozenbin * stimPara.refreshrate/stimPara.Nblinks;

lrRsq            = NaN(Ncells,1);  lrCCnorm    = NaN(Ncells,1);
modelRsq         = NaN(Ncells,1);  modelCCnorm = NaN(Ncells,1);

lrpredictions    = NaN(Ncells, stimPara.FrozenFrames-Nt+1, 'single');
modelpredictions = NaN(Ncells, stimPara.FrozenFrames-Nt+1, 'single');

% filter with temporal components and predict
for icell = 1:Ncells
    
    cellrates = squeeze(trialrates(icell,:,:));
    % low rank
    fgens_lr = conv(filteredstim_lr(icell,:), flip(temporalComponents(icell,:)),'valid');
    preds_lr = rf.getPredictionFromBinnedNonlinearity(fgens_lr,...
        nlncentslr(icell, :), nlnvalslr(icell, :));
    
    % model
    fgens_model = conv(filteredstim_model(icell,:), flip(modeltcomps(icell,:)),'valid');
    preds_model = rf.getPredictionFromBinnedNonlinearity(fgens_model,...
        nlncentsmodel(icell, :), nlnvalsmodel(icell, :));
    
    lrCCnorm(icell)    = rf.calc_CCnorm(cellrates', preds_lr');
    lrRsq(icell)       = rf.rsquare(mean(cellrates, 2), preds_lr');
    modelCCnorm(icell) = rf.calc_CCnorm(cellrates', preds_model');
    modelRsq(icell)    = rf.rsquare(mean(cellrates, 2), preds_model');
    modelpredictions(icell,:) = preds_model;
    lrpredictions   (icell,:) = preds_lr;
    
end

res.lrpredictions = lrpredictions;          res.modelpredictions    = modelpredictions;
res.lrRsq         = lrRsq;                  res.modelRsq            = modelRsq;
res.lrCCnorm      = lrCCnorm;               res.modelCCnorm         = modelCCnorm;
res.stimPara      = stimPara;
fprintf('Done! Took %2.2f s...\n', toc);

end


function res = nlpreddichromatic(res, runningbin, frozenbin, stimPara)

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;
Ncols = numel(stimPara.seed);
Ncells = size(runningbin,1);
Ntrials = size(res.trialRates,3);

%--------------------------------------------------------------------------
fprintf('Calculating generator signals... '); tic;

% first separate and normalize the tempoal component
tcs = res.temporalComponents;
tcs = [squeeze(tcs(:,:,1)), squeeze(tcs(:,:,2))];
tcsnorm = tcs./sqrt(sum(tcs.^2,2));
tccol1 = tcsnorm(:,1:Nt);
tccol2 = tcsnorm(:,Nt+1:end);

% then the spatial component
spc = res.spatialComponents;
spc = [squeeze(spc(:,:,:,1)), squeeze(spc(:,:,:,2))];
lrspacepredict  = reshape(flip(spc, 2), Ncells, Ny*Nx*2);
lrspacepredict  = lrspacepredict./sqrt(sum(lrspacepredict.^2,2));
% this is filp of the orignial scs
lrspacepredict = reshape(lrspacepredict,Ncells,Ny*Ncols,Nx);
lrspacepredictcol1 = lrspacepredict(:,1:Ny,:);
lrspacepredictcol1 = reshape(lrspacepredictcol1, Ncells, Ny*Nx);
lrspacepredictcol2 = lrspacepredict(:,Ny+1:end,:);
lrspacepredictcol2 = reshape(lrspacepredictcol2, Ncells, Ny*Nx);

% calculate low-rank generator signals
lrgeneratorscol1 = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    lrspacepredictcol1, tccol1, stimPara.seed(1));
lrgeneratorscol2 = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    lrspacepredictcol2, tccol2, stimPara.seed(2));

% first separate and normalize the tempoal component
tcs = res.modeltcomps;
tcs = [squeeze(tcs(:,:,1)), squeeze(tcs(:,:,2))];
tcsnorm = tcs./sqrt(sum(tcs.^2,2));
modeltccol1 = tcsnorm(:,1:Nt);
modeltccol2 = tcsnorm(:,Nt+1:end);

% then the spatial component
spc = res.modelscomps;
spc = [squeeze(spc(:,:,:,1)), squeeze(spc(:,:,:,2))];
modelspredict  = reshape(flip(spc, 2), Ncells, Ny*Nx*2);
modelspredict  = modelspredict./sqrt(sum(modelspredict.^2,2));
% this is filp of the orignial scs
modelspredict = reshape(modelspredict,Ncells,Ny*Ncols,Nx);
modelspredictcol1 = modelspredict(:,1:Ny,:);
modelspredictcol1 = reshape(modelspredictcol1, Ncells, Ny*Nx);
modelspredictcol2 = modelspredict(:,Ny+1:end,:);
modelspredictcol2 = reshape(modelspredictcol2, Ncells, Ny*Nx);

% calculate model generator signals
modelgeneratorscol1 = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    modelspredictcol1, modeltccol1, stimPara.seed(1));
[modelgeneratorscol2, spikes] = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    modelspredictcol2, modeltccol2, stimPara.seed(2));


fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
fprintf('Extracting nonlinearities...\n'); tic;

nltypes = {'combined', 'conditional','nl2d','marginal'};

for ii = 1: numel(nltypes)
    
    [ nl_lowrank.([nltypes{ii},'_nly']), nl_lowrank.([nltypes{ii},'_nlx'])] = rfroutines.getDichromaticNonlinearities(nltypes{ii},...
        lrgeneratorscol1, lrgeneratorscol2, spikes, stimPara.nonlinBinN, stimPara.Nblinks/stimPara.refreshrate);

    [ nl_model.([nltypes{ii},'_nly']), nl_model.([nltypes{ii},'_nlx'])] = rfroutines.getDichromaticNonlinearities(nltypes{ii},...
        modelgeneratorscol1, modelgeneratorscol2, spikes, stimPara.nonlinBinN, stimPara.Nblinks/stimPara.refreshrate);
    fprintf('calculating %s nonlinearities! Took %.2f s...\n', nltypes{ii}, toc);

end

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
% for low-rank nonlinearities
lrpredictions = rfroutines.getDichromaticPredictions(frozenbin, nl_lowrank, lrspacepredictcol1, ...
    lrspacepredictcol2, tccol1, tccol2, stimPara);
% for model nonlinearities
modelpredictions = rfroutines.getDichromaticPredictions(frozenbin, nl_model, modelspredictcol1, ...
    modelspredictcol2, modeltccol1, modeltccol2, stimPara);

% setting up the output
res.nl_lowrank          =   nl_lowrank;
res.nl_model            =   nl_model;
res.lrpredictions       =   lrpredictions;
res.modelpredictions    =   modelpredictions;


end


function res = nlpredtrichromatic(res, runningbin, frozenbin, stimPara)


end

