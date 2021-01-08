

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

[stimPara, ftimes, spikes, clusters, savingpath] = getRFparameters(datapath);

[res, spiketimes, runningbin, frozenbin] = calcRFsta(ftimes, spikes, clusters, stimPara);

res = rfTemporalSpatialComps(res, spiketimes, clusters, stimPara);

rfdata = rfroutines.rfNonlinearitiesPredictions(res, runningbin, frozenbin, stimPara);
rfdata.savingpath = savingpath;
rfdata.clusters = clusters;
rfdata.para = stimPara;

if not(isempty(savingpath))
    fprintf('Saving data... '); tic;
    saveName = [filesep, num2str(stimPara.expnumber,'%02d'), '-checkerflicker_analysis_for_experiment_on_',stimPara.date,'.mat'];
    save([savingpath,filesep,'rf_data',saveName], '-v7.3', '-struct', 'rfdata');
    fprintf('Done! Took %2.2f s...\n', toc);
end

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


function [stimPara, ft, spiketimes, clus, savingpath] = getRFparameters(dp, varargin)

expData = loadRawData(dp,'frozennoise','RF_Analysis','rf_data','excludename','1600x1200'); % only load data, no folder making

stimPara = expData.stimPara;
%stimPara = stimdata.stimPara;
stimPara = checkstimParaArgs(stimPara, 'stimulus', 'frozennoise');
stimPara = checkstimParaArgs(stimPara, 'refreshrate', 60);
stimPara = checkstimParaArgs(stimPara, 'screen', expData.screen);
stimPara = checkstimParaArgs(stimPara, 'fs', double(expData.samplingrate));
stimPara = checkstimParaArgs(stimPara, 'pulseRate', 2);

% for auto-correlogram
stimPara = checkstimParaArgs(stimPara, 'dtcorr', 5e-4);         % 5e-4;
stimPara = checkstimParaArgs(stimPara, 'Ncorr', 60e-3 / 5e-4); % %60e-3/dtcorr; % old values
stimPara = checkstimParaArgs(stimPara, 'normACG', true);

stimPara.lightprojection = expData.lightprojection;
if strcmpi( expData.lightprojection, 'oled')
    pixsize = 7.5e-6;
else % now for lightcrafter, add option for patch setups later
    if isfield(expData,'recordingsetup')
        recsetup = expData.recordingsetup;
    else
        recsetup = questdlg('Select the recording setup?','Recording setup','Aragorn','Bilbo','Aragorn');
    end
    switch lower(recsetup)
        case 'aragorn'
            pixsize = 7.2e-6;
            stimPara.recordingsetup = 'Aragorn';
        case 'bilbo'
            pixsize = 7.6e-6;
            stimPara.recordingsetup = 'Aragorn';
        otherwise
            error('Which setup you used for recording, cabron, how do get the pixel size?');
    end
end
stimPara = checkstimParaArgs(stimPara, 'pixelsize', pixsize);

% this is to get optimal number of bins for sta and nonlinearity
upsampres = floor(1e3/stimPara.refreshrate*stimPara.nblinks);      % ==> up-sampling ratio, change this for different resolution.
timebehindspikeinSec = 680;   %  700 ms behind each spike
resbyFrame = round(floor(1e3/stimPara.refreshrate*stimPara.nblinks)) / upsampres;
nonlinBinN = floor(timebehindspikeinSec/(1e3/stimPara.refreshrate*stimPara.nblinks)) *resbyFrame;
if (nonlinBinN/resbyFrame) < 15, nonlinBinN = 15 * resbyFrame; end % minimum bins are 15

stimPara = checkstimParaArgs(stimPara, 'filterWindow', timebehindspikeinSec/ 1e3);
stimPara = checkstimParaArgs(stimPara, 'nsigma', 2);
stimPara = checkstimParaArgs(stimPara, 'nonlinBinN', nonlinBinN);

if ~isfield(stimPara,'color')
    colormode = 'monochromatic';
else
    if ~stimPara.color
        colormode = 'monochromatic';
    elseif ~stimPara.usered && stimPara.usegreen && stimPara.useblue
        colormode = 'dichromatic';
    elseif stimPara.usered && stimPara.usegreen && stimPara.useblue
        colormode = 'trichromatic';
    end
end
stimPara.colormode = colormode;

switch colormode
    case 'monochromatic'
        stimPara.seed = stimPara.seedrunningnoise;
        stimPara.secondseed = stimPara.seedfrozennoise;
        stimPara.usedcolor = false;
        stimPara.colororder = {'blackwhite'};
        stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise'});
        
    case 'dichromatic'
        stimPara.seed = [stimPara.seedrunninggreen, stimPara.seedrunningblue];
        stimPara.secondseed = [ stimPara.seedfrozengreen, stimPara.seedfrozenblue];
        stimPara.meanintensity = [stimPara.greenmeanintensity, stimPara.bluemeanintensity];
        stimPara.contrast = [stimPara.greenContrast, stimPara.blueContrast];
        stimPara.usedcolor = [stimPara.usered, stimPara.usegreen, stimPara.useblue];
        stimPara.colororder = {'green','blue/UV'};
        stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise','seedrunningred','seedrunninggreen',...
            'seedrunningblue','seedfrozenred','seedfrozengreen','seedfrozenblue','redmeanintensity',...
            'greenmeanintensity','bluemeanintensity','redContrast','greenContrast','blueContrast',...
            'usered','usegreen', 'useblue'});
        
    case 'trichromatic'
        stimPara.seed = [stimPara.seedrunningred, stimPara.seedrunninggreen, stimPara.seedrunningblue];
        stimPara.secondseed = [stimPara.seedfrozenred, stimPara.seedfrozengreen, stimPara.seedfrozenblue];
        stimPara.meanintensity = [stimPara.redmeanintensity,  stimPara.greenmeanintensity, stimPara.bluemeanintensity];
        stimPara.contrast = [stimPara.redContrast, stimPara.greenContrast, stimPara.blueContrast];
        stimPara.usedcolor = [stimPara.usered, stimPara.usegreen, stimPara.useblue];
        stimPara.colororder = {'red','green','blue/UV'};
        stimPara = rmfield(stimPara,{'seedrunningnoise','seedfrozennoise','seedrunningred','seedrunninggreen',...
            'seedrunningblue','seedfrozenred','seedfrozengreen','seedfrozenblue','redmeanintensity',...
            'greenmeanintensity','bluemeanintensity','redContrast','greenContrast','blueContrast',...
            'usered','usegreen', 'useblue'});
end

chsize = [num2str(stimPara.stixelheight),'x',num2str(stimPara.stixelwidth)];    % checking stixel size
if stimPara.nblinks == 1,   bltxt  = 'blink';   else,       bltxt = 'blinks';   end
savingpath = [dp,filesep,'Data Analysis',filesep,num2str(stimPara.expnumber,'%02d'),'-FrozenNoise_',...
    colormode,'_',chsize,'_',num2str(stimPara.nblinks),bltxt];

%%%folder making
if ~exist(savingpath,'dir'), mkdir(savingpath); end
if ~exist([savingpath,'/rf_data'],'dir')
    mkdir ([savingpath,'/rf_data']);
end


if ~iscolumn(expData.ftimes), expData.ftimes = expData.ftimes'; end

if size(expData.ftimes,2) ~=2
    if stimPara.nblinks == 1 && max(size(expData.ftimes)) > (max(size(expData.ftimesoff))*2-10)
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

if isfield(stimPara,'nblinks'), stimPara.Nblinks = stimPara.nblinks; stimPara = rmfield(stimPara,'nblinks'); end
stimPara.Nx = ceil(stimPara.screen(1)/stimPara.stixelwidth);
stimPara.Ny = ceil(stimPara.screen(2)/stimPara.stixelheight);

stimPara.date = expData.date;
if isfield(expData,'sortinginfo')
    stimPara.sortinfo = expData.sortinginfo;
end

spiketimes = expData.spiketimes;
clus = expData.clusters;

end


function para = checkstimParaArgs(para, argname, defval)
if not(isfield(para, argname))
    para.(argname) = defval;
end
end


function [res, spiketimes, runningbin, frozenbin] = calcRFsta(ftimes, spikes, clus, stimPara, varargin)

disp(['Starting ' stimPara.stimulus ' analysis for stimulus ' num2str(stimPara.expnumber) '...'])

Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;
Ncols = numel(stimPara.seed);

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
res.trialRates = frozenbin;

allReliableRsq = rfroutines.imageTrialRsq( permute(frozenbin,[1 3 2]) ); %think of removing variable part
res.allReliableRsq = allReliableRsq;

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

stas = zeros(Ncells, Ny, Nx, Nt, Ncols,'single');
for ii = 1:numel(stimPara.seed)
    staAll = rfroutines.calculateBlockSTAbwGPU(runningbin, Nt, Nx*Ny, stimPara.seed(ii), stimPara.contrast(ii));
    staAll = reshape(staAll, Ncells, Ny, Nx,Nt);
    staAll = flip(staAll,2); %flipping because of C++
    stas(:,:,:,:,ii) = staAll;
end

res.stas = squeeze(stas);
fprintf('Done! \n');
%--------------------------------------------------------------------------
%define coordinates
timeVec=(-(Nt-1/2):1:-1/2)*stimPara.Nblinks/stimPara.refreshrate; %in seconds
spaceVecX = stimPara.lmargin +0.5 + stimPara.stixelwidth*(0:Nx-1)+stimPara.stixelwidth/2;
spaceVecY = stimPara.bmargin +0.5 + stimPara.stixelheight*(0:Ny-1)+stimPara.stixelheight/2;
res.timeVec = timeVec;
res.spaceVecX = spaceVecX;
res.spaceVecY = spaceVecY;
%--------------------------------------------------------------------------
end


function res = rfTemporalSpatialComps(res, spiketimes, clus, stimPara, varargin)


Nt = ceil(stimPara.filterWindow*stimPara.refreshrate/stimPara.Nblinks);
Ny = stimPara.Ny; Nx = stimPara.Nx;
Ncols = numel(stimPara.seed);
Ncells = size(clus,1);
nelipspoints = 100;
timeVec=(-(Nt-1/2):1:-1/2)*stimPara.Nblinks/stimPara.refreshrate; %in seconds
rfac   = 4.5 * 1.4826;

stixelsForFit = ceil(50/stimPara.stixelwidth);

gaussParams        = NaN(Ncells, 6, Ncols);
spatialComponents  = zeros(Ncells, Ny, Nx, Ncols,'single');
modelscomps        = zeros(Ncells, Ny, Nx, Ncols,'single');
temporalComponents = zeros(Ncells, Nt, Ncols, 'single');
modeltcomps        = zeros(Ncells, Nt, Ncols,'single');
autoCorrs          = NaN(Ncells, stimPara.Ncorr, 'single');
rfdiameters        = NaN(Ncells, Ncols);
ellipseareas       = NaN(Ncells, Ncols);
contourareas       = NaN(Ncells, Ncols);
contourpoints      = cell(Ncells, Ncols);
ellipsepoints      = NaN(Ncells, 2, nelipspoints, Ncols, 'single');
rfmodelparams      = NaN(Ncells, 12, Ncols);
allrangex          = cell(Ncells, Ncols);
allrangey          = cell(Ncells, Ncols);
allmoran           = NaN(Ncells, Ncols);
allmads            = NaN(Ncells, Ncols);
iuse               = NaN(Ncells, Ncols);
rfcolresp          = NaN(Ncells,1);

cellspktimes = accumarray(spiketimes(:,2),spiketimes(:,1), [Ncells, 1], @(x) {x});
dpx = stimPara.pixelsize;


for icol = 1:Ncols
    staAll = squeeze(res.stas(:,:,:,:,icol));
    %Nstfit = 2*stixelsForFit + 1;
    %allspred = zeros(Ncells, 1);
    allmads(:,icol) = mad(staAll(:,:),1,2);
    %rsta = staAll(:,:)./allmads;
    iuse(:,icol) = sum((abs(staAll(:,:)) - rfac*allmads(:,icol))>0, 2);
end
%contfac = 0.14;
%dtcorr  = stimPara.dtcorr;  %5e-4;
%Ncorr   = stimPara.Ncorr;   %60e-3/dtcorr; % old values

disp('Beginning cell by cell analysis...'); tic;
msg = [];
for icell = 1:Ncells
    
    %zoomsta = zeros((stixelsForFit*2+1).^2, Nt, Ncols);
    sigtcpix = cell(Ncols, 1);
    % we loop twice, first to find which color is better and the cut the region
    % around the rf from that receptive field.
    
    for icol = 1: Ncols
        if iuse(icell, icol) == 0, continue; end
        csta = double(squeeze(res.stas(icell,:,:,:,icol)));
        %======================================================================
        % select ROI after blurring
        smsta = rfroutines.smoothSTA(squeeze(res.stas(icell,:,:,:,icol)), 0.5);
        [~, imax] = max(abs(smsta(:)));
        [y0, x0, ~] = ind2sub(size(smsta), imax);
        
        rangeX = x0+(-stixelsForFit:stixelsForFit); rangeX = rangeX(rangeX>0 & rangeX<=Nx);
        rangeY = y0+(-stixelsForFit:stixelsForFit); rangeY = rangeY(rangeY>0 & rangeY<=Ny);
        zoomsta = reshape(csta(rangeY, rangeX, :), numel(rangeY)*numel(rangeX), Nt);
        
        allrangex{icell,icol} = rangeX;
        allrangey{icell,icol} = rangeY;
        
        %======================================================================
        % find significant pixels in the zoomed region
        [bpx, ~] = find(abs(zoomsta) > rfac*mad(csta(:),1));
        if isempty(bpx), continue, end
        sigtcpix{icol} = bpx;
        %======================================================================
        % extract temporal and spatial components
        
        tempcomp = mean(zoomsta(bpx,:),1)';
        spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));
        allmoran(icell,icol) = rfroutines.moransI(spcomp, size(spcomp,1), size(spcomp,2));
        
    end
    %[s,d,v] = svd(zoomsta,'econ');
    
    [~, bestcolrf] = max(allmoran(icell,:)); % select the best color response
    rfcolresp(icell) = bestcolrf;
    
    for icol = 1: Ncols
        
        if iuse(icell, icol) == 0, continue; end
        csta = double(squeeze(res.stas(icell,:,:,:,icol)));
        
        rangeX  = allrangex{icell, bestcolrf};
        rangeY  = allrangey{icell, bestcolrf};
        zoomsta = reshape(csta(rangeY, rangeX, :), numel(rangeY)*numel(rangeX), Nt);
        
        % correct it back to best color
        allrangex{icell,icol} = rangeX;
        allrangey{icell,icol} = rangeY;
        
        %======================================================================
        % find significant pixels in the zoomed region
        %[bpx, ~] = find(abs(zoomsta) > rfac*mad(csta(:),1));
        bpx = sigtcpix{bestcolrf};
        if isempty(bpx), continue, end
        %======================================================================
        % extract temporal and spatial components
        
        tempcomp = mean(zoomsta(bpx,:),1)';
        spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));
        
        
        spatialComponents(icell, rangeY,  rangeX, icol) = spcomp;
        temporalComponents(icell, :, icol) = tempcomp;
        
        
        % get simple ellipse fit and contour
        
        %tempcomp = mean(zoomsta(bpx,:),1)';
        %spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));
        
        %allmoran(icell) = rfroutines.moransI(spcomp, size(spcomp,1), size(spcomp,2));
        
        %spatialComponents(icell, rangeY,  rangeX) = spcomp;
        %temporalComponents(icell, :) = tempcomp;
        
        % get simple ellipse fit and contour
        
        [contpts, contarea, centgaussparams] = rfroutines.getRfContourPts(...
            res.spaceVecX(rangeX), res.spaceVecY(rangeY), spcomp);
        
        
        
        %     c = getEllipseFromParams(centgaussparams, 2);
        %     clf;
        %     imagesc(spaceVecX(rangeX),spaceVecY(rangeY), spcomp, [-1 1]*max(abs(spcomp(:))));
        %     hold on; colormap(redblue)
        %     plot(contpts(1,:), contpts(2,:), '-k', c(1,:), c(2,:), 'g')
        %
        %     areac = contarea*(dpx*1e3)^2;
        %     effd  = 2*sqrt(areac/pi);
        %     title(sprintf('diam: %d', round(effd*1e3)))
        
        gaussParams(icell, :, icol) = centgaussparams;
        contourareas(icell, icol)   = contarea*(dpx*1e3)^2;
        rfdiam = rfroutines.getRFDiam(rfroutines.getGaussFromParams(centgaussparams), stimPara.nsigma, dpx);
        ellipseareas(icell,icol) = pi * (rfdiam*1e3/2)^2;
        contourpoints{icell, icol} = contpts;
        ellipsepoints(icell,:,:,icol) = rfroutines.getEllipseFromParams(centgaussparams, stimPara.nsigma, nelipspoints);
        
        
        %==========================================================================
        % fit DoG+time RF model
        tempGuess  = rfroutines.fitTemporalComponent(timeVec, tempcomp);
        spaceGuess = rfroutines.dogreceptfield2(res.spaceVecX(rangeX), res.spaceVecY(rangeY), spcomp);
        
        fullGuess = [tempGuess spaceGuess(1:6) spaceGuess(8)/spaceGuess(7)];
        stafit    = permute(reshape(zoomsta,numel(rangeY),numel(rangeX), Nt), [3 1 2]);
        modelprms = rfroutines.fitParametricSTA2(timeVec, res.spaceVecX(rangeX), res.spaceVecY(rangeY), stafit, fullGuess);
        rfmodelparams(icell, :, icol) = modelprms;
        
        modeltcomps(icell, :, icol) = rfroutines.templowpass(modelprms(1:5), timeVec);
        
        spParams = [modelprms(6:end-1) 1 modelprms(end)];
        [X, Y] = meshgrid(res.spaceVecX(rangeX),res.spaceVecY(rangeY));
        Z = rfroutines.dogmatrixfun(spParams,{X(:), Y(:)});
        modelscomps(icell, rangeY,  rangeX, icol) = reshape(Z, numel(rangeY), numel(rangeX));
        
        
        %==========================================================================
        % get model values
        %tvec = templowpass(modelprms(1:5), timeVec);
        
        cgauss = rfroutines.getGaussFromParams(modelprms(6:end));
        rfdiameters(icell, icol) = rfroutines.getRFDiam(cgauss, stimPara.nsigma, dpx);
        
    end
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


bigestcont = max(max(cellfun(@(x) (size(x,2)),contourpoints)));
contpts = nan(Ncells,2,bigestcont+1, Ncols, 'single');
surrIdx = nan(Ncells, Ncols);
for icol = 1: Ncols
    cp = contourpoints(:,icol);
    for ii = 1:Ncells
        if not(isempty(cp{ii}))
            contpts(ii,:,1:size(cp{ii},2),icol) = cp{ii};
        end
    end
    
    sigmas = linspace(0,8,1e3);
    %central equation
    activations = (1-exp(-sigmas.^2/2))-...
        rfmodelparams(:,end,icol).*(1-exp(-sigmas.^2/2./rfmodelparams(:,11,icol).^2));
    activations(activations<0) = 0;
    surrIdx(:,icol) = 1-activations(:,end)./max(activations,[],2);
    
end
% autocorr output
if stimPara.normACG
    autoCorrs = autoCorrs./sum(autoCorrs,2);
end
autoCorrsLag   = linspace(0, stimPara.Ncorr * stimPara.dtcorr *1e3, stimPara.Ncorr); % xaxis of autocorr
%--------------------------------------------------------------------------
res.spatialComponents  = spatialComponents;         res.modelscomps     = modelscomps;
res.temporalComponents = temporalComponents;        res.modeltcomps     = modeltcomps;
res.autoCorrelations   = autoCorrs;                 res.autoCorrLag     = autoCorrsLag;
res.surroundIdx     = surrIdx;
res.sigmaActivation    = activations;               res.sigmaVals       = sigmas;
res.allrangex          = allrangex;                 res.allrangey       = allrangey;
res.contourareas       = contourareas;              res.ellipseareas    = ellipseareas;
res.rfdiameters        = rfdiameters;               res.contourpoints   = contpts;
res.gaussparams        = gaussParams;               res.rfmodelparams   = rfmodelparams;
res.ellipsepoints      = ellipsepoints;             res.allmoran        = allmoran;
res.rfcolresp          = rfcolresp;

fn = fieldnames(res);
for ii = 1:numel(fn)
    if ndims (res.(fn{ii}) ) >= 3 %#ok
        res.(fn{ii}) = squeeze(res.(fn{ii}));
    end
end

end

