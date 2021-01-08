

function res = getDichromaticPredictions(frozenbin, nl,spcompCol1, spcompCol2, tcompCol1, tcompCol2, stimPara, varargin)

fprintf('Calculating predictions and performances... ');
tic;

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;
Ncols = numel(stimPara.seed);
Ncells = size(frozenbin,1);

% first color
frozenstimCol1 = rfroutines.ran1bool(stimPara.secondseed(1), stimPara.FrozenFrames*Nx*Ny);
frozenstimCol1 = reshape(frozenstimCol1, Ny*Nx, stimPara.FrozenFrames);
frozenstimCol1 = 2 * single(frozenstimCol1) - 1; % transform to contrast

% second color
frozenstimCol2 = rfroutines.ran1bool(stimPara.secondseed(2), stimPara.FrozenFrames*Nx*Ny);
frozenstimCol2 = reshape(frozenstimCol2, Ny*Nx, stimPara.FrozenFrames);
frozenstimCol2 = 2 * single(frozenstimCol2) - 1; % transform to contrast

% filter with spatial components
filteredstim_col1 = spcompCol1  * frozenstimCol1;
filteredstim_col2 = spcompCol2  * frozenstimCol2;

% getting trial rates
trialrates = frozenbin * stimPara.refreshrate/stimPara.Nblinks;

% pre-allocating memory
[nlcombined_CCnorm, nlcombined_Rsq] = deal(nan(Ncells,1));
nlcombined_predictions = nan(Ncells, stimPara.FrozenFrames-Nt+1, 'single');
[nlcond_CCnorm, nlcond_Rsq] = deal(nan(Ncells,Ncols));
nlcond_predictions = nan(Ncells, stimPara.FrozenFrames-Nt+1, Ncols, 'single');
[nlmarginal_CCnorm, nlmarginal_Rsq] = deal(nan(Ncells,Ncols));
nlmarginal_predictions = nan(Ncells, stimPara.FrozenFrames-Nt+1, Ncols, 'single');

% filter with temporal components and predict
for icell = 1:Ncells
    
    cellrates = squeeze(trialrates(icell,:,:));
    
    fgens_Col1 = conv(filteredstim_col1(icell,:), flip(tcompCol1(icell,:)),'valid');
    fgens_Col2 = conv(filteredstim_col2(icell,:), flip(tcompCol2(icell,:)),'valid');
    
    if all(isnan(fgens_Col1)) || all(isnan(fgens_Col2)) ,continue;  end
    
    
    % first combined nonlinearities
    fgens = fgens_Col1 + fgens_Col2;
    preds_combi = rfroutines.getPredictionFromBinnedNonlinearity(fgens,...
        nl.combined_nlx(icell, :), nl.combined_nly(icell, :));
    
    nlcombined_CCnorm(icell) = rfroutines.calc_CCnorm(cellrates', preds_combi');
    nlcombined_Rsq(icell)    = rfroutines.rsquare(mean(cellrates, 2), preds_combi');
    nlcombined_predictions(icell,:) = preds_combi;
    
    
    % second conditional nonlinearities
    preds_cond1 = rfroutines.getPredictionFromBinnedNonlinearity(fgens_Col1,...
        squeeze(nl.conditional_nlx(icell, 1, :)), squeeze(nl.conditional_nly(icell, 1, :)));
    
    preds_cond2 = rfroutines.getPredictionFromBinnedNonlinearity(fgens_Col2,...
        squeeze(nl.conditional_nlx(icell, 2, :)), squeeze(nl.conditional_nly(icell, 2, :)));
    
    nlcond_CCnorm(icell,1) = rfroutines.calc_CCnorm(cellrates', preds_cond1');
    nlcond_CCnorm(icell,2) = rfroutines.calc_CCnorm(cellrates', preds_cond2');
    nlcond_Rsq(icell,1)    = rfroutines.rsquare(mean(cellrates, 2), preds_cond1');
    nlcond_Rsq(icell,2)    = rfroutines.rsquare(mean(cellrates, 2), preds_cond2');
    nlcond_predictions(icell,:,1) = preds_cond1;
    nlcond_predictions(icell,:,2) = preds_cond2;
    
    
    % third marginal nonlinearities
    preds_marg1 = rfroutines.getPredictionFromBinnedNonlinearity(fgens_Col1,...
        squeeze(nl.marginal_nlx(icell, 1, :)), squeeze(nl.marginal_nly(icell, 1, :)));
    
    preds_marg2 = rfroutines.getPredictionFromBinnedNonlinearity(fgens_Col2,...
        squeeze(nl.marginal_nlx(icell, 2, :)), squeeze(nl.marginal_nly(icell, 2, :)));
    
    nlmarginal_CCnorm(icell,1) = rfroutines.calc_CCnorm(cellrates', preds_marg1');
    nlmarginal_CCnorm(icell,2) = rfroutines.calc_CCnorm(cellrates', preds_marg2');
    nlmarginal_Rsq(icell,1)    = rfroutines.rsquare(mean(cellrates, 2), preds_marg1');
    nlmarginal_Rsq(icell,2)    = rfroutines.rsquare(mean(cellrates, 2), preds_marg2');
    nlmarginal_predictions(icell,:,1) = preds_marg1;
    nlmarginal_predictions(icell,:,2) = preds_marg2;
    
end

res.combined_predictions       = nlcombined_predictions;
res.combined_CCnorm            = nlcombined_CCnorm;
res.combined_Rsq               = nlcombined_Rsq;

res.conditional_predictions    = nlcond_predictions;
res.conditional_CCnorm         = nlcond_CCnorm;
res.conditional_Rsq            = nlcond_Rsq;

res.marginal_predictions       = nlmarginal_predictions;
res.marginal_CCnorm            = nlmarginal_CCnorm;
res.marginal_Rsq               = nlmarginal_Rsq;

fprintf('Done! Took %2.2f s...\n', toc);

end