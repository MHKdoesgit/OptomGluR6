

function varargout = calc_FrozenRF(datapath, varargin)
%
%%% calc_FrozenRF %%%
%
%
% This function calculate spatio-temporal receptive field of ganglion cells
% based on thier reponse to runnig-frozen checker flicker stimulus
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   RFData : a cell structure containing stimulus parameters,2D STAs, norm of STA,
%           fitted 2D STAs, svd of STAs, svd of zoom STAs, best spatial and
%           temporal compenents, and fitted spatial component with 2D Gaussian
%           function with parameters of fit.
%   Plot : This function will also plot the svd, fitted spatial and
%          temporal components and receptive field area of ganglion cell.
%
% written by Mohammad based on calc_RF function on 25.10.2018.

totaltime =tic;
if (nargin < 1), datapath = uigetdir(); end    % for no input opens a uigetdir
%if (nargin > 1), cell2analyze = varargin{1}; else,  cell2analyze = 0; end
% get all the STAs and other parameters

[stimPara, ftimes, spikes, clusters, savingpath] = rf_analysis_parameters(datapath);

[res, spiketimes, runningbin, frozenbin] = calc_rf_sta(ftimes, spikes, clusters, stimPara);
rfdata = rf_tempcomp_spcomp_nonlin(res, spiketimes, runningbin, frozenbin, clusters, stimPara, savingpath);

% plotting all the receptive fields together
% plotFrozennoise(RFdata, para, clusters, savingpath)

disp(' And BOOM!!! Goes all the Receptive Fields ');
disp(seconds2human (toc(totaltime)));
sound(struct2array(load('chirp.mat','y')))
varargout{1} = rfdata;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [stimPara, ft, spiketimes, clus, savingpath] = rf_analysis_parameters(dp, varargin)

expData = loadRawData(dp,'frozennoise','RF_Analysis','rf_data','excludename','1600x1200'); % only load data, no folder making

stimPara = expData.stimPara;
expinfo  = expData.info;
scinfo   = expinfo.screen;
%stimPara = stimdata.stimPara;
stimPara = check_stimpara_args(stimPara, 'stimulus', 'frozennoise');
stimPara = check_stimpara_args(stimPara, 'refreshrate', scinfo.refreshrate);
stimPara = check_stimpara_args(stimPara, 'screen', scinfo.resolution);
stimPara = check_stimpara_args(stimPara, 'fs', double(expinfo.samplingrate));
stimPara = check_stimpara_args(stimPara, 'pulseRate', 2);

% for auto-correlogram
stimPara = check_stimpara_args(stimPara, 'dtcorr', 5e-4);         % 5e-4;
stimPara = check_stimpara_args(stimPara, 'Ncorr', 60e-3 / 5e-4); % %60e-3/dtcorr; % old values
stimPara = check_stimpara_args(stimPara, 'normACG', true);

stimPara = check_stimpara_args(stimPara, 'lightprojection', scinfo.type);
% pixsize  = scinfo.pixelsize * 1e-6;
% if strcmpi( expData.lightprojection, 'oled')
%     pixsize = 7.5e-6;
% else % now for lightcrafter, add option for patch setups later
%     if isfield(expData,'recordingsetup')
%         recsetup = expData.recordingsetup;
%     else
%         recsetup = questdlg('Select the recording setup?','Recording setup','Aragorn','Bilbo','Aragorn');
%     end
%     switch lower(recsetup)
%         case 'aragorn'
%             pixsize = 7.2e-6;
%             stimPara.recordingsetup = 'Aragorn';
%         case 'bilbo'
%             pixsize = 7.6e-6;
%             stimPara.recordingsetup = 'Aragorn';
%         otherwise
%             error('Which setup you used for recording, cabron, how do get the pixel size?');
%     end
% end
stimPara = check_stimpara_args(stimPara, 'pixelsize', scinfo.pixelsize * 1e-6);
if isfield(stimPara,'nblinks'), stimPara.Nblinks = stimPara.nblinks; stimPara = rmfield(stimPara,'nblinks'); end


% this is to get optimal number of bins for sta and nonlinearity
upsampres = floor(1e3/stimPara.refreshrate*stimPara.Nblinks);      % ==> up-sampling ratio, change this for different resolution.
timebehindspikeinSec = 600;   % ms behind each spike
resbyFrame = round(floor(1e3/stimPara.refreshrate*stimPara.Nblinks)) / upsampres;
nonlinBinN = floor(timebehindspikeinSec/(1e3/stimPara.refreshrate*stimPara.Nblinks)) *resbyFrame;
if (nonlinBinN/resbyFrame) < 15, nonlinBinN = 15 * resbyFrame; end % minimum bins are 15
% nonlinBinN = ceil(nonlinBinN ./ 5)*5;

stimPara = check_stimpara_args(stimPara, 'filterWindow', timebehindspikeinSec/ 1e3);
stimPara = check_stimpara_args(stimPara, 'nsigma', 2);
stimPara = check_stimpara_args(stimPara, 'nonlinBinN', nonlinBinN);

if isfield(stimPara,'seedrunningnoise'), stimPara.seed = stimPara.seedrunningnoise; stimPara = rmfield(stimPara,'seedrunningnoise'); end
if isfield(stimPara,'seedfrozennoise'), stimPara.secondseed = stimPara.seedfrozennoise; stimPara = rmfield(stimPara,'seedfrozennoise'); end
stimPara = check_stimpara_args(stimPara, 'colororder', 'blackwhite');
% % remove margins if they are zero
% if stimPara.lmargin == 0, stimPara = rmfield(stimPara, 'lmargin'); end
% if stimPara.rmargin == 0, stimPara = rmfield(stimPara, 'rmargin'); end
% if stimPara.tmargin == 0, stimPara = rmfield(stimPara, 'tmargin'); end
% if stimPara.bmargin == 0, stimPara = rmfield(stimPara, 'bmargin'); end



% if ~isfield(stimPara,'color')
%     colormode = 'monochromatic';
% else
%     if ~stimPara.color
%         colormode = 'monochromatic';
%     elseif ~stimPara.usered && stimPara.usegreen && stimPara.useblue
%         colormode = 'dichromatic';
%     elseif stimPara.usered && stimPara.usegreen && stimPara.useblue
%         colormode = 'trichromatic';
%     end
% end
% stimPara.colormode = colormode;
%
% switch colormode
%     case 'monochromatic'
%         stimPara.seed = stimPara.seedrunningnoise;
%         stimPara.secondseed = stimPara.seedfrozennoise;
%         stimPara.usedcolor = false;
%         stimPara.colororder = {'blackwhite'};
%         stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise'});
%
%     case 'dichromatic'
%         stimPara.seed = [stimPara.seedrunninggreen, stimPara.seedrunningblue];
%         stimPara.secondseed = [ stimPara.seedfrozengreen, stimPara.seedfrozenblue];
%         stimPara.meanintensity = [stimPara.greenmeanintensity, stimPara.bluemeanintensity];
%         stimPara.contrast = [stimPara.greenContrast, stimPara.blueContrast];
%         stimPara.usedcolor = [stimPara.usered, stimPara.usegreen, stimPara.useblue];
%         stimPara.colororder = {'green','blue/UV'};
%         stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise','seedrunningred','seedrunninggreen',...
%             'seedrunningblue','seedfrozenred','seedfrozengreen','seedfrozenblue','redmeanintensity',...
%             'greenmeanintensity','bluemeanintensity','redContrast','greenContrast','blueContrast',...
%             'usered','usegreen', 'useblue'});
%
%     case 'trichromatic'
%         stimPara.seed = [stimPara.seedrunningred, stimPara.seedrunninggreen, stimPara.seedrunningblue];
%         stimPara.secondseed = [stimPara.seedfrozenred, stimPara.seedfrozengreen, stimPara.seedfrozenblue];
%         stimPara.meanintensity = [stimPara.redmeanintensity,  stimPara.greenmeanintensity, stimPara.bluemeanintensity];
%         stimPara.contrast = [stimPara.redContrast, stimPara.greenContrast, stimPara.blueContrast];
%         stimPara.usedcolor = [stimPara.usered, stimPara.usegreen, stimPara.useblue];
%         stimPara.colororder = {'red','green','blue/UV'};
%         stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise','seedrunningred','seedrunninggreen',...
%             'seedrunningblue','seedfrozenred','seedfrozengreen','seedfrozenblue','redmeanintensity',...
%             'greenmeanintensity','bluemeanintensity','redContrast','greenContrast','blueContrast',...
%             'usered','usegreen', 'useblue'});
% end

chsize = [num2str(stimPara.stixelheight),'x',num2str(stimPara.stixelwidth)];    % checking stixel size
if stimPara.Nblinks == 1,   bltxt  = 'blink';   else,       bltxt = 'blinks';   end
savingpath = [dp,filesep,'Data Analysis',filesep,num2str(stimPara.expnumber,'%02d'),'-FrozenNoise_',...
    chsize,'_',num2str(stimPara.Nblinks),bltxt];

%%%folder making
if ~exist(savingpath,'dir'), mkdir(savingpath); end
% if ~exist([savingpath,'/rf_data'],'dir')
%     mkdir ([savingpath,'/rf_data']);
%end


if ~iscolumn(expData.ftimes), expData.ftimes = expData.ftimes'; end

if size(expData.ftimes,2) ~=2
    if stimPara.Nblinks == 1 && max(size(expData.ftimes)) > (max(size(expData.ftimesoff))*2-10)
        ft = [expData.ftimes(1:2:numel(expData.ftimes)),expData.ftimes(2:2:numel(expData.ftimes))];
    else
        if isrow(expData.ftimes) && isrow(expData.ftimesoff)
            ft = [expData.ftimes(1:end-1)',expData.ftimesoff'];
        else
            ft = [expData.ftimes,expData.ftimesoff];
        end
    end
else
    ft = expData.ftimes;
end



stimPara.Nx = ceil(stimPara.screen(1)/stimPara.stixelwidth);
stimPara.Ny = ceil(stimPara.screen(2)/stimPara.stixelheight);

stimPara.date = expData.date;
if isfield(expData,'sortinginfo')
    stimPara.sortinfo = expData.sortinginfo;
end

% if isfield(expData,'info')
%     stimPara.expinfo = expData.info;
% end

spiketimes = expData.spiketimes;
clus = expData.clusters;

end

%--------------------------------------------------------------------------------------------------%

function para = check_stimpara_args(para, argname, defval)
if not(isfield(para, argname))
    para.(argname) = defval;
end
end

%--------------------------------------------------------------------------------------------------%

function [res, spiketimes, runningbin, frozenbin] = calc_rf_sta(ftimes, spikes, clus, stimPara, varargin)

disp(['Starting ' stimPara.stimulus ' analysis for stimulus ' num2str(stimPara.expnumber) '...'])

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;

Ncells = size(clus,1);
if iscell(spikes) && size(spikes,2)~=2
    spiketimes = rfroutines.spikeCell2Mat(spikes,stimPara.fs);
    spiketimes = spiketimes(:,[2 1]); % the dimension of spiketimes should be flipped to match KS order
    ft = round(ftimes * stimPara.fs);
else
    spiketimes = spikes;
    ft = ftimes;
end
spikesbin  = rfroutines.blinkBinnerKS(spiketimes, Ncells, ft(:,1), ft(:,2), stimPara.Nblinks, stimPara.pulseRate);

%spikesbin = stimdata.spikesbin;
Nframes   = size(spikesbin,2);

runningFrames = stimPara.RunningFrames;
trialFrames   = runningFrames+stimPara.FrozenFrames;

Ntrials     = floor(Nframes/trialFrames);
totalFrames = Ntrials*trialFrames;

totalbin   = reshape(spikesbin(:,1:totalFrames), Ncells, trialFrames, Ntrials);
runningbin = totalbin(:,1:runningFrames,:);

%this part adds the remaining spikes for STA calculation
rembin     = spikesbin(:,totalFrames+1:end);
if size(rembin,2)>runningFrames; rembin=rembin(:,1:runningFrames); end
runningbin = cat(3,runningbin,zeros(Ncells,runningFrames));
runningbin(:,1:size(rembin,2),end) = rembin;
%--------------------------------------------------------------------------
% Prediction part
frozenbin = totalbin(:, runningFrames+Nt:end,:);
res.frozenTrialRates = frozenbin;

allReliableRsq = rfroutines.imageTrialRsq( permute(frozenbin,[1 3 2]) ); %think of removing variable part
res.frozenRsq = allReliableRsq;

frozenRates = mean(frozenbin,3)*stimPara.refreshrate/stimPara.Nblinks;
res.frozenRates = frozenRates;

frozenTimeVec = (0:(stimPara.FrozenFrames-1))*stimPara.Nblinks/stimPara.refreshrate; %in seconds
frozenTimeVec = frozenTimeVec+stimPara.Nblinks/stimPara.refreshrate/2;
res.frozenTimeVec = frozenTimeVec;
%--------------------------------------------------------------------------
disp('Generating STAs for the running part...');
%seeduse = stimPara.seed;
% if isfield(stimPara, 'color')
%     if stimPara.color
%         seeduse = stimPara.seed/2;
%     end
% end

%sta = zeros(Ncells, Ny, Nx, Nt,'single');
staAll = rfroutines.calculateBlockSTAbwGPU(runningbin, Nt, Nx*Ny, stimPara.seed, stimPara.contrast);
staAll = reshape(staAll, Ncells, Ny, Nx,Nt);
res.sta  = flip(staAll,2); %flipping because of C++
%    sta(:,:,:,:) = staAll;

%res.sta = squeeze(sta);
fprintf('Done! \n');
%--------------------------------------------------------------------------
%define coordinates
timeVec=(-(Nt-1/2):1:-1/2)*stimPara.Nblinks/stimPara.refreshrate; %in seconds
spaceVecX = stimPara.lmargin +0.5 + stimPara.stixelwidth*(0:Nx-1)+stimPara.stixelwidth/2;
spaceVecY = stimPara.bmargin +0.5 + stimPara.stixelheight*(0:Ny-1)+stimPara.stixelheight/2;
res.tcTimeVec = timeVec;
res.spaceVecX = spaceVecX;
res.spaceVecY = spaceVecY;

end

%--------------------------------------------------------------------------------------------------%

function res = rf_tempcomp_spcomp_nonlin(res, spiketimes, runningbin, frozenbin, clus, stimPara, savingpath, varargin)

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny;
Nx = stimPara.Nx;
Ncells = size(clus,1);
Ntrials = size(res.frozenTrialRates,3);
nelipspoints = 100;
timeVec=(-(Nt-1/2):1:-1/2)*stimPara.Nblinks/stimPara.refreshrate; %in seconds
rfac   = 4.5 * 1.4826;
stixelsForFit = ceil(50/stimPara.stixelwidth);
spaceVecX = res.spaceVecX;      spaceVecY = res.spaceVecY;
%--------------------------------------------------------------------------

allmads = mad(res.sta(:,:),1,2);
iuse = sum((abs(res.sta(:,:)) - rfac*allmads)>0, 2);

%contfac = 0.14;
%dtcorr  = stimPara.dtcorr;  %5e-4;
%Ncorr   = stimPara.Ncorr;   %60e-3/dtcorr; % old values

gaussParams        = NaN(Ncells, 6);
spatialComponents  = zeros(Ncells, Ny, Nx,'single');
modelscomps        = zeros(Ncells, Ny, Nx,'single');
temporalComponents = zeros(Ncells, Nt,'single');
modeltcomps        = zeros(Ncells, Nt,'single');
autoCorrs          = NaN(Ncells, stimPara.Ncorr, 'single');
rfdiameters        = NaN(Ncells, 1);
ellipseareas       = NaN(Ncells, 1);
contourareas       = NaN(Ncells, 1);
contourpoints      = cell(Ncells, 1);
ellipsepoints      = NaN(Ncells, 2, nelipspoints);
rfmodelparams      = NaN(Ncells, 12);
allrangex          = cell(Ncells, 1);
allrangey          = cell(Ncells, 1);
allmoran           = NaN(Ncells, 1);
% get spike times
cellspktimes = accumarray(spiketimes(:,2),spiketimes(:,1), [Ncells, 1], @(x) {x});

disp('Beginning cell by cell analysis...'); tic;
msg = [];
for icell = 1:Ncells
    
    if iuse(icell) == 0, continue; end
    csta = double(squeeze(res.sta(icell,:,:,:)));
    %======================================================================
    % select ROI after blurring
    smsta = rfroutines.smoothSTA(squeeze(res.sta(icell,:,:,:)), 0.5);
    [~, imax] = max(abs(smsta(:)));
    [y0, x0, ~] = ind2sub(size(smsta), imax);
    
    rangeX = x0+(-stixelsForFit:stixelsForFit); rangeX = rangeX(rangeX>0 & rangeX<=Nx);
    rangeY = y0+(-stixelsForFit:stixelsForFit); rangeY = rangeY(rangeY>0 & rangeY<=Ny);
    zoomsta = reshape(csta(rangeY, rangeX, :), numel(rangeY)*numel(rangeX), Nt);
    
    allrangex{icell} = rangeX; allrangey{icell} = rangeY;
    %======================================================================
    % find significant pixels in the zoomed region
    [bpx, ~] = find(abs(zoomsta) > rfac*mad(csta(:),1));
    if isempty(bpx), continue, end
    %======================================================================
    % extract temporal and spatial components
    
    tempcomp = mean(zoomsta(bpx,:),1)';
    spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));
    
    %[s,d,v] = svd(zoomsta,'econ');
    
    spatialComponents(icell, rangeY,  rangeX) = spcomp;
    temporalComponents(icell, :) = tempcomp;
    
    % get simple ellipse fit and contour    
    allmoran(icell) = rfroutines.moransI(spcomp, size(spcomp,1), size(spcomp,2));
    
    % get simple ellipse fit and contour    
    [contpts, contarea, centgaussparams] = rfroutines.getRfContourPts(...
        spaceVecX(rangeX),spaceVecY(rangeY), spcomp);
    
    gaussParams(icell, :) = centgaussparams;
    contourareas(icell)   = contarea*(stimPara.pixelsize*1e3)^2;
    rfdiam = rfroutines.getRFDiam(rfroutines.getGaussFromParams(centgaussparams), stimPara.nsigma, stimPara.pixelsize);
    ellipseareas(icell) = pi * (rfdiam*1e3/2)^2;
    contourpoints{icell} = contpts;
    ellipsepoints(icell,:,:) = rfroutines.getEllipseFromParams(centgaussparams, stimPara.nsigma, nelipspoints);
    %==========================================================================
    % fit DoG+time RF model
    tempGuess  = rfroutines.fitTemporalComponent(timeVec, tempcomp);
    spaceGuess = rfroutines.dogreceptfield2(spaceVecX(rangeX), spaceVecY(rangeY), spcomp);
    
    fullGuess = [tempGuess spaceGuess(1:6) spaceGuess(8)/spaceGuess(7)];
    stafit    = permute(reshape(zoomsta,numel(rangeY),numel(rangeX), Nt), [3 1 2]);
    modelprms = rfroutines.fitParametricSTA2(timeVec, spaceVecX(rangeX), spaceVecY(rangeY), stafit, fullGuess);
    rfmodelparams(icell, :) = modelprms;
    
    modeltcomps(icell, :) = rfroutines.templowpass(modelprms(1:5), timeVec);
    
    spParams = [modelprms(6:end-1) 1 modelprms(end)];
    [X, Y] = meshgrid(spaceVecX(rangeX),spaceVecY(rangeY));
    Z = rfroutines.dogmatrixfun(spParams,{X(:), Y(:)});
    modelscomps(icell, rangeY,  rangeX) = reshape(Z, numel(rangeY), numel(rangeX));
    %==========================================================================
    % get model values
    %tvec = templowpass(modelprms(1:5), timeVec);
    
    cgauss = rfroutines.getGaussFromParams(modelprms(6:end));
    rfdiameters(icell) = rfroutines.getRFDiam(cgauss, stimPara.nsigma, stimPara.pixelsize);
    %==========================================================================
    % get acg
    spks = cellspktimes{icell}/stimPara.fs;
    K = rfroutines.ccg(spks, spks, stimPara.Ncorr, stimPara.dtcorr);
    K(stimPara.Ncorr+1) = 0;
    autoCorrs(icell, :) = single(K((stimPara.Ncorr+1):end-1));
    %==========================================================================
    if mod(icell, 20) == 0 || icell == Ncells
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', icell,Ncells,toc);
        fprintf(msg);
    end
end

% to change the format of contourpoints from cell to mat for faster plotting
cp = contourpoints;
bigestcont = max(cellfun(@(x) (size(x,2)),cp));
contpts = nan(Ncells,2,bigestcont+1,'single');
for ii = 1:Ncells
    if not(isempty(cp{ii}))
        contpts(ii,:,1:size(cp{ii},2)) = cp{ii};
    end
end

sigmas = linspace(0,8,1e3);
%central equation
activations = (1-exp(-sigmas.^2/2))-...
    rfmodelparams(:,end).*(1-exp(-sigmas.^2/2./rfmodelparams(:,11).^2));
activations(activations<0) = 0;
surrIdx = 1-activations(:,end)./max(activations,[],2);
% autocorr output
if stimPara.normACG
    autoCorrs = autoCorrs./sum(autoCorrs,2);
end
autoCorrsLag   = linspace(0, stimPara.Ncorr * stimPara.dtcorr *1e3, stimPara.Ncorr); % xaxis of autocorr
%--------------------------------------------------------------------------
res.spatialComps        = spatialComponents;         res.modelspcomps       = modelscomps;
res.temporalComps       = temporalComponents;        res.modeltcomps        = modeltcomps;
res.autoCorrelations    = autoCorrs;                 res.autoCorrLag        = autoCorrsLag;
res.surroundIndex       = surrIdx;                   res.surrActivation     = activations;
res.surrSigmaVals       = sigmas;
res.rangex              = allrangex;                 res.rangey             = allrangey;
res.contourareas        = contourareas;              res.ellipseareas       = ellipseareas;
res.rfdiameters         = rfdiameters;               res.contourpoints      = contpts;
res.gaussparams         = gaussParams;               res.rfmodelparams      = rfmodelparams;
res.ellipsepoints       = ellipsepoints;             res.moransI            = allmoran;

%--------------------------------------------------------------------------
fprintf('Calculating generator signals... '); tic;

temporalComponents = temporalComponents./sqrt(sum(temporalComponents.^2,2));
lrspacepredict  = reshape(flip(spatialComponents, 2), Ncells, Ny*Nx);
lrspacepredict  = lrspacepredict./sqrt(sum(lrspacepredict.^2,2));

% calculate low-rank generator signals
[lrgenerators, ~] = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
    lrspacepredict, temporalComponents, stimPara.seed);

modeltcomps  = modeltcomps./sqrt(sum(modeltcomps.^2,2));
modelspredict  = reshape(flip(modelscomps, 2), Ncells, Ny*Nx);
modelspredict  = modelspredict./sqrt(sum(modelspredict.^2,2));

% calculate model generator signals
[modelgenerators, spikes] = rfroutines.calculateLowRankGeneratorsbw(runningbin(:,:,1:Ntrials),...
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
    [staVals, staCents,~] = rfroutines.getNonlinearity(cellgenslr, spikes(icell,:),...
        stimPara.nonlinBinN,stimPara.Nblinks/stimPara.refreshrate);
    nlncentslr(icell,:) = staCents;
    nlnvalslr(icell, :) = staVals;
    
    cellgensmodel = modelgenerators(icell,:)';
    if isnan(sum(cellgensmodel)); continue; end
    [staVals,staCents,~] = rfroutines.getNonlinearity(cellgensmodel, spikes(icell,:),...
        stimPara.nonlinBinN,stimPara.Nblinks/stimPara.refreshrate);
    nlncentsmodel(icell,:) = staCents;
    nlnvalsmodel(icell,:) = staVals;
end

res.nlx         = nlncentslr;               res.nly             = nlnvalslr;
res.nlxmodel    = nlncentsmodel;            res.nlymodel        = nlnvalsmodel;

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
fprintf('Calculating predictions and performances... '); tic;

%Generate stimulus
[frozenstimulus, ~] = rfroutines.ran1bool(stimPara.secondseed, stimPara.FrozenFrames*Nx*Ny);
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
    preds_lr = rfroutines.getPredictionFromBinnedNonlinearity(fgens_lr,...
        nlncentslr(icell, :), nlnvalslr(icell, :));
    
    % model
    fgens_model = conv(filteredstim_model(icell,:), flip(modeltcomps(icell,:)),'valid');
    preds_model = rfroutines.getPredictionFromBinnedNonlinearity(fgens_model,...
        nlncentsmodel(icell, :), nlnvalsmodel(icell, :));
    
    lrCCnorm(icell)    = rfroutines.calc_CCnorm(cellrates', preds_lr');
    lrRsq(icell)       = rfroutines.rsquare(mean(cellrates, 2), preds_lr');
    modelCCnorm(icell) = rfroutines.calc_CCnorm(cellrates', preds_model');
    modelRsq(icell)    = rfroutines.rsquare(mean(cellrates, 2), preds_model');
    modelpredictions(icell,:) = preds_model;
    lrpredictions   (icell,:) = preds_lr;
    
end

res.predictions     = lrpredictions;          res.modelpredictions      = modelpredictions;
res.predRsq         = lrRsq;                  res.modelpredRsq          = modelRsq;
res.predCCnorm      = lrCCnorm;               res.modelpredCCnorm       = modelCCnorm;

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
res.savingpath = savingpath;
res.clusters = clus;
res.para = stimPara;
rfdata = res;
if not(isempty(savingpath))
    fprintf('Saving data... '); tic;
    saveName = [filesep, num2str(stimPara.expnumber,'%02d'), '-checkerflicker_analysis_for_experiment_on_',stimPara.date,'.mat'];
    save([savingpath,saveName], '-v7.3', '-struct', 'rfdata');
    fprintf('Done! Took %2.2f s...\n', toc);
end
%--------------------------------------------------------------------------


end
