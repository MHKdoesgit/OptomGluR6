

function varargout = Analyze_FullFieldFlicker(datapath, varargin)
%
%%% Analyze_FrozenNoise %%%
%
%
% This function analyzes the FrozenNoise stimulus. It first calculates the
% rasters and psth for the frozen parts of the stimulus and then calculated
% the STA for the running part. After that it will uses the sta and
% nonlinearities of the runing noise to predict the forzen noise parts. The
% output is the predication of the frozen part along with the stas and
% nonlinearities of the runing parts.
%
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   fndataall : a cell structure containing rasters values, psth values and
%             stas and nonlinearities.
%   Plot : This function will also plot the rasters and psth and the stas
%          and nonlinearities together.
%
% written by Mohammad, 20.07.2016
% updated to completely new version with faster calculation for stas and
% nonlinearities on 28.05.2018.
% added the option for plottig of monochromatic frozen running noise on
% 12.11.2018.


totaltime =tic;
if (nargin < 1),  datapath = uigetdir();       end

% get the stimulus separated into running and frozen parts and get the
% frozen part analysis
tic;
[spikbins,frozenspikes,para,clusters,savingpath] = fff_parameters(datapath,0);
ffdata = calcSTAandstimEnsembles(spikbins,frozenspikes,clusters, para, savingpath);
toc;
% if strcmpi(para.stimMode, 'monochromatic')
%     plotFrozenNoiseMonochromatic(fndataall, para, clusters, savingpath)
% elseif strcmpi(para.stimMode, 'dichromatic')
%     plotFrozenNoiseDichromatic(fndataall, ciff, para, clusters, savingpath);
% end
sound(struct2array(load('gong.mat','y')));
disp(seconds2human (toc(totaltime)));
varargout{1} = ffdata;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [spkbin, frozenspikes, stimPara, clus, savingpath] = fff_parameters(dp,varargin)

[expData, savingpath] = loadRawData(dp,'fullfieldfrozennoise','FullFieldFlicker_Analysis');
% loading clusters and frametimings
clus = expData.clusters;
ft = expData.ftimes;

stimPara = expData.stimPara;
expinfo  = expData.info;
scinfo   = expinfo.screen;
stimPara.refreshrate = scinfo.refreshrate;
stimPara.screen = scinfo.resolution;
stimPara.fs = double(expinfo.samplingrate);
stimPara.pulseRate = 2;
stimPara.lightprojection = scinfo.type;
stimPara.dtcorr = 5e-4;
stimPara.Ncorr = 60e-3 / 5e-4;
stimPara.normACG = true;
if isfield(stimPara,'nblinks'), stimPara.Nblinks = stimPara.nblinks; stimPara = rmfield(stimPara,'nblinks'); end

stimPara.significancecheck = false; % STC significance check (very slow!)


stimPara.upsamplescale = 1;
% if isempty((regexpi(dp,{'color','marmoset'})))
%     stimPara.stimMode = 'trichromatic';
% else
%     stimPara.stimMode = 'dichromatic';
% end % checking color mode
% if ~any(isfield(expData.stimPara,{'color','Color'})), stimPara.stimMode = 'monochromatic'; end
% if stimPara.color == 0, stimPara.stimMode = 'monochromatic'; end

% this is to get optimal number of bins for sta and nonlinearity
upsampres = floor(1e3/stimPara.refreshrate*stimPara.Nblinks);      % ==> up-sampling ratio, change this for different resolution.
timebehindspikeinSec = 600;   % ms behind each spike
resbyFrame = round(floor(1e3/stimPara.refreshrate*stimPara.Nblinks)) / upsampres;
nonlinBinN = floor(timebehindspikeinSec/(1e3/stimPara.refreshrate*stimPara.Nblinks)) *resbyFrame;
if (nonlinBinN/resbyFrame) < 15, nonlinBinN = 15 * resbyFrame; end % minimum bins are 15
stimPara.filterWindow = timebehindspikeinSec/ 1e3;
stimPara.nonlinBinN = nonlinBinN;
if isfield(stimPara,'seedrunningnoise'), stimPara.seed = stimPara.seedrunningnoise; stimPara = rmfield(stimPara,'seedrunningnoise'); end
if isfield(stimPara,'seedfrozennoise'), stimPara.secondseed = stimPara.seedfrozennoise; stimPara = rmfield(stimPara,'seedfrozennoise'); end

% remove margins if they are zero
if stimPara.lmargin == 0, stimPara = rmfield(stimPara, 'lmargin'); end
if stimPara.rmargin == 0, stimPara = rmfield(stimPara, 'rmargin'); end
if stimPara.tmargin == 0, stimPara = rmfield(stimPara, 'tmargin'); end
if stimPara.bmargin == 0, stimPara = rmfield(stimPara, 'bmargin'); end

%stimPara.numSTA = (40 / stimPara.nblinks) * stimPara.resolution;
%stimPara.numNL = stimPara.numSTA;
stimPara.nRepeats =  floor(length(ft) / (stimPara.RunningFrames + stimPara.FrozenFrames));
%stimPara.screen = expData.screen;
% stimPara = rmfield(stimPara,{'originalname','expnumber','stimulus','meanintensity','contrast',...
%     'lmargin','rmargin','bmargin','tmargin'});
% if strcmpi(stimPara.stimMode,'dichromatic')
%     stimPara = rmfield(stimPara,{'redmeanintensity','greenmeanintensity','bluemeanintensity'});
% end

ftblinks = round(mean(diff(ft))*stimPara.refreshrate);
ft = ft(1:stimPara.Nblinks/ftblinks:end);
stimPara.nRepeats = ceil(stimPara.nRepeats/(stimPara.Nblinks/ftblinks));
numfullft = stimPara.nRepeats *  (stimPara.RunningFrames + stimPara.FrozenFrames);
missingft = numfullft-length(ft);
ft = [ft, max(ft) + interp1(1:length(ft),ft,1:missingft,'linear')];
ft = ft(1:stimPara.nRepeats *  (stimPara.RunningFrames + stimPara.FrozenFrames));
if stimPara.upsamplescale > 1
    ftHR = linspace(ft(1),ft(end),length(ft) * stimPara.upsamplescale);
else
    ftHR = ft;
end

spkbin = nan(length(expData.spiketimes),length(ftHR));

for ii= 1:length(expData.spiketimes)
    if ~isempty(expData.spiketimes{ii})
        spkbin(ii,:) = (histc(expData.spiketimes{ii},ftHR)); %#ok
        % / ((1/60*para.nblinks)/para.resolution);
    end
end

stimPara.date = expData.date;
if isfield(expData,'sortinginfo')
    stimPara.sortinfo = expData.sortinginfo;
end

frozenspikes = getFrozenNoiseSpikesPSTH(ft, expData.spiketimes,stimPara);
end

%--------------------------------------------------------------------------------------------------%

function frnspikes = getFrozenNoiseSpikesPSTH(ft,spk,para,varargin)

%ft = ft(1:((para.RunningFrames + para.FrozenFrames)*para.nRepeats));
ft2D = reshape(ft,[],para.nRepeats);
frnft = ft2D(para.RunningFrames+1:(para.RunningFrames + para.FrozenFrames),:);
%frzduration = para.FrozenFrames * (1/(60/para.nblinks));
%frzpsth = zeros(size(spk,2),para.FrozenFrames-para.numSTA+1);
frnspikes = cell(size(spk,2),1);
for ii = 1:length(spk)
    thisspk = spk{ii}';
    allfrnspk = cell(1,para.nRepeats);
    for jj = 1:para.nRepeats
        allfrnspk{jj} = thisspk(and(thisspk >= frnft(1,jj),thisspk <= frnft(para.FrozenFrames,jj))) - frnft(1,jj);
    end
    frnspikes{ii} = CelltoMatUE(allfrnspk);
    % thispsth = mean(histc(frnspikes{ii},linspace(0,frzduration,para.FrozenFrames+1),1),2);
    % frzpsth(ii,:) = thispsth(para.numSTA:para.FrozenFrames) / (1/(60/para.nblinks));
    %clear thisspk allfrnspk;
end

end

%--------------------------------------------------------------------------------------------------%

function res = calcSTAandstimEnsembles(spikesbin,frozenspikes, clus, para, savingpath)

rnseed = para.seed;
frseed = para.secondseed;
% res = cell(numel(rnseed),1);

Ncells=size(clus,1);
trialFrames = para.RunningFrames+para.FrozenFrames;
spikesbin=reshape(spikesbin(:,1:para.nRepeats*trialFrames),Ncells,trialFrames,[]);
toHzratio = para.Nblinks/para.refreshrate; % normalize the response to Hz
% signumiter = 100; % number iteration for STC significance testing


%--------------------------------------------------------------------------
disp('Analyzing the frozen parts...');



frozenstimulus = gasdev(frseed, para.FrozenFrames)';
frozenHankel= hankel(1:para.nonlinBinN,para.nonlinBinN:para.FrozenFrames);
res.frozenEn = frozenstimulus(frozenHankel);

frozenbin = spikesbin(:,para.RunningFrames+1:end,:);
res.trialrates = frozenbin(:,para.nonlinBinN:end,:) * para.refreshrate/para.Nblinks;
res.frozenRates = mean(frozenbin(:,para.nonlinBinN:end,:),3)*para.refreshrate/para.Nblinks;
res.frozentvec=(para.nonlinBinN:(para.FrozenFrames))*toHzratio+toHzratio/2;
res.frozenspikes = frozenspikes;
%--------------------------------------------------------------------------
disp('Analyzing the running parts');

runningbin=spikesbin(:,1:para.RunningFrames,:);

seed = rnseed;%para.seedrunningblue;%stimPara.seed;
rnstimEn=zeros(para.nRepeats,para.nonlinBinN, para.RunningFrames-para.nonlinBinN+1);
for itrial = 1:para.nRepeats
    [stimulus,seed] =   gasdev(seed, para.RunningFrames);
    stimulus        =   stimulus';
    hankelMat       =   hankel(1:para.nonlinBinN,para.nonlinBinN:para.RunningFrames);
    rnstimEn(itrial,:,:)    =   stimulus(hankelMat);
end
runningbin=reshape(permute(runningbin(:,para.nonlinBinN:end,:),[1 3 2]),Ncells,[]);
rnstimEn=reshape(permute(rnstimEn, [2 1 3]), para.nonlinBinN,[]);

res.sta = bsxfun(@rdivide,runningbin*rnstimEn',sum(runningbin,2));
%--------------------------------------------------------------------------

res.staNorm = bsxfun(@times,res.sta,1./sqrt(sum(res.sta.^2,2))); %STAs are normalized!!!

res.runnoiseEn = rnstimEn;
res.runnoisespkbins = runningbin;

onr=abs(max(res.staNorm,[],2)); offr=abs(min(res.staNorm,[],2));
res.biphasicIndex=(onr-offr)./(onr+offr);
%res.biphasicIndex=biphasicIndex;

res.timeVec=(-(para.nonlinBinN-1):1:0)*toHzratio; %in seconds
%res.timeVec=timeVec;

res.runningGenerators = res.staNorm * rnstimEn;
res.frozenGenerators = res.staNorm * res.frozenEn;

modelsta = zeros(size(res.sta));
for ii = 1:Ncells
    tempGuess  = rfroutines.fitTemporalComponent(res.timeVec, res.staNorm(ii,:)');
    modelsta(ii,:) = rfroutines.templowpass(tempGuess, res.timeVec);
end
res.modelsta = modelsta;
res.modelrunGenerators = res.modelsta * rnstimEn;
res.modelfrozGenerators = res.modelsta * res.frozenEn;
res = calcnonlinprediction(res,para);
% STC for each color (without significance test! )
res.stc = getffSTC(res,para);


%clearvars frozenstimulus frozenHankel frozenbin runningbin seed itrial rnstimEn stimulus hankelMat onr offr;
res.savingpath = savingpath;
res.clusters = clus;
res.para = para;
ffdata = res;
if not(isempty(savingpath))
    fprintf('Saving data... '); tic;
    saveName = [filesep, num2str(para.expnumber,'%02d'), '-fullfieldflicker_analysis_for_experiment_on_',para.date,'.mat'];
    save([savingpath,saveName], '-v7.3', '-struct', 'ffdata');
    fprintf('Done! Took %2.2f s...\n', toc);
end


% fndataall = cell2mat(fndataall);

end

%--------------------------------------------------------------------------------------------------%

function res = calcnonlinprediction(res,para)

Ncell       =   size(res.runnoisespkbins,1);
%nlx         =   nan(Ncell, para.nonlinBinN);
%nly         =   nan(Ncell, para.nonlinBinN);
[predictions, modelpreds]     =   deal(nan(size(res.frozenRates)));
%predstat    =   struct([]);
toHzratio   =   para.Nblinks/para.refreshrate; % normalize the response to Hz
[nlx, nly, nlxmodel, nlymodel] = deal(nan(Ncell, para.nonlinBinN));
[Rsqr, pearsonCoeff, explVar, CCnorm, Rsqrmodel, pearsonCoeffmodel, explVarmodel, CCnormmodel] = deal(nan(Ncell,1));

%allEigenvalues=NaN(cellN, windowN);
%allSTCNCvecs=NaN(cellN, windowN);
%allSTCNCcents=NaN(cellN, options.nonlinBinN);
%allSTCNCvals=NaN(cellN, options.nonlinBinN);
%predictionsSTCNC=NaN(size(frozenRates));
%rsqSTCNC=NaN(cellN,1);

for ii=1:Ncell
    
    if all(res.frozenRates(ii,:)==0)
        continue;
    end
    cellspikes = res.runnoisespkbins(ii,:);
    %stagenerators=res.runningGenerators(ii,:);
    [nly(ii,:),nlx(ii,:),~] = rfroutines.getNonlinearity(res.runningGenerators(ii,:)', cellspikes', para.nonlinBinN,toHzratio);
    %nlx(ii,:)=stacents;
    %nly(ii,:)=stavals;
    
    [nlymodel(ii,:),nlxmodel(ii,:)] = rfroutines.getNonlinearity(res.modelrunGenerators(ii,:)', cellspikes', para.nonlinBinN,toHzratio);
    
    
    %--------------------------------------------------------------------------
    %     STC=getSTC(stimEn', cellspikes,0);
    %     if sum(isinf(STC(:)))|| sum(isnan(STC(:))); continue; end;
    %
    %     [stcVecs,stcEigs]=eig(STC);
    %     [stcEigs,sortedInds]=sort(diag(stcEigs));
    %     stcVecs=stcVecs(:,sortedInds);
    %     allEigenvalues(ii,:)=stcEigs;
    %     stcncVec=stcVecs(:,end);
    %     stcncgenerators=stcncVec'*stimEn;
    %     [stcncvals,stcnccents,~]=getNonlinearity(stcncgenerators', cellspikes',...
    %         options.nonlinBinN,stimPara.Nblinks/experiment.projector.refreshrate);
    %     if stcncvals(1)>stcncvals(end)
    %         stcncvals=flipud(stcncvals);
    %         stcnccents=-fliplr(stcnccents);
    %         stcncVec=-stcncVec;
    %     end
    %     allSTCNCvecs(ii,:)=stcncVec;
    %     allSTCNCcents(ii,:)=stcnccents;
    %     allSTCNCvals(ii,:)=stcncvals;
    %--------------------------------------------------------------------------
    if sum(nly(ii,:)>0)>0
        predictions(ii,:) = rfroutines.getPredictionFromBinnedNonlinearity(res.frozenGenerators(ii,:), nlx(ii,:), nly(ii,:));
        modelpreds(ii,:) = rfroutines.getPredictionFromBinnedNonlinearity(res.frozenGenerators(ii,:), nlxmodel(ii,:),nlymodel(ii,:));
        
        %predictions(ii,:)=predsSTA;
        
        [Rsqr(ii),pearsonCoeff(ii),explVar(ii)] = calcRsqPearsonCoeff(res.frozenRates(ii,:),predictions(ii,:));
        
        [Rsqrmodel(ii),pearsonCoeffmodel(ii),explVarmodel(ii)] = calcRsqPearsonCoeff(res.frozenRates(ii,:),modelpreds(ii,:));
        
        cellrates = squeeze(res.trialrates(ii,:,:));
        
        CCnorm(ii)    = rfroutines.calc_CCnorm(cellrates',predictions(ii,:)');
        CCnormmodel(ii)    = rfroutines.calc_CCnorm(cellrates',modelpreds(ii,:)');
        %%lrRsq(icell)       = rfroutines.rsquare(mean(cellrates, 2), predsSTA');
        
        
        
        %        rsqSTA(ii)=rsquare(frozenRates(ii,:),predsSTA);
        %
        %         frozengensSTCNC = stcncVec'*frozenEn;
        %         predsSTCNC = getPredictionFromBinnedNonlinearity(frozengensSTCNC, stcnccents, stcncvals);
        %         predictionsSTCNC(ii,:)=predsSTCNC;
        %         rsqSTCNC(ii)=rsquare(frozenRates(ii,:),predsSTCNC);
    end
end
%%
res.nlx = nlx;
res.nly = nly;
res.predictions = predictions;
res.nlxmodel = nlxmodel;
res.nlymodel = nlymodel;
res.modelpredictions = modelpreds;
res.stat.predRsq = Rsqr;
res.stat.pearsonCoeff = pearsonCoeff;
res.stat.explVar = explVar;
res.stat.predCCnorm = CCnorm;
res.stat.modelpredRsq = Rsqrmodel;
res.stat.modelpearsonCoeff = pearsonCoeffmodel;
res.stat.modelexplVar = explVarmodel;
res.stat.modelpredCCnorm = CCnormmodel;

end

%--------------------------------------------------------------------------------------------------%

function stcout = getffSTC(res,para, varargin)

if nargin > 2, numiter = varargin{1}; else, numiter = 100; end
%binnedspks = binnedspks';
Ncells = size(res.runnoisespkbins,1);
nlbin = para.nonlinBinN;
toHzratio = para.Nblinks/para.refreshrate; % normalize the response to Hz
covMat = getSTC(res.runnoiseEn, res.runnoisespkbins');
peakreg = ceil(prctile(1:nlbin,80));
[~, maxpos] = (max(abs(res.staNorm(:,peakreg:end)),[],2));
peakstasign = floor(diag(res.staNorm(:,(peakreg + maxpos -1))))*2+1;
stcout = cell(Ncells,1);
fprintf('Checking significant eigenvalues and associated eigenvectors... '); tic;

for ii = 1:Ncells
    if all(all(isnan(squeeze (covMat(:,:,ii)))))    % this is to avoid cases with no spikes
        [stcout{ii}.stceigs, stcout{ii}.sigvecs, stcout{ii}.sigvecsnorm, stcout{ii}.nlx, stcout{ii}.nly] = deal(zeros(nlbin,6));
        stcout{ii}.stcvecs = zeros(nlbin,nlbin);
        stcout{ii}.sigeigs = zeros(6,1);
        stcout{ii}.confint = zeros(2,nlbin);
        stcout{ii}.eigminidx = 0;
        stcout{ii}.eigmaxidx = nlbin;
        stcout{ii}.significancecheck = false;
        stcout{ii} = orderfields(stcout{ii},[1 6:11 2 3 5 4]);
        warning(['Ayayaya, we dont got no spikes here for channel ',num2str(ii)]);
        continue;
    end
    
    [stcvecs,stceigs] = eig(squeeze (covMat(:,:,ii)));
    [stcout{ii}.stceigs, sortidx] = sort(diag(stceigs));
    stcout{ii}.stcvecs = stcvecs(:,sortidx);
    
    if not(isreal(stcout{ii}.stcvecs )), stcout{ii}.stcvecs = real(stcout{ii}.stcvecs); end
    if not(isreal(stcout{ii}.stceigs )), stcout{ii}.stceigs = real(stcout{ii}.stceigs); end
    
    if para.significancecheck    % check for significance of eigen values by comparing it to random set (slower)
        [sigvecs, stcout{ii}.sigeigs, stcout{ii}.confint, stcout{ii}.eigminidx, stcout{ii}.eigmaxidx] = ...
            findSigEigs(res.runnoiseEn', res.runnoisespkbins(ii,:), numiter, 0.01, 0);
        stcout{ii}.significancecheck = true;
    end
    
    if ~para.significancecheck || isempty(sigvecs)  % instead of checking for significance just get the first and last 3 eigen values
        sigvecs = [stcout{ii}.stcvecs(:,1:3),stcout{ii}.stcvecs(:,end-2:end)];
        stcout{ii}.sigeigs = [stcout{ii}.stceigs(1:3);stcout{ii}.stceigs(end-2:end)];
        stcout{ii}.confint = zeros(2,nlbin);
        stcout{ii}.eigminidx = 3;
        stcout{ii}.eigmaxidx = nlbin-3;
        stcout{ii}.significancecheck = false;
        %disp(['no significant shit was found for channel ',num2str(ii)]);
    end
    
    % correct the sign of the eigen vector to match the sign of the STA
    for jj = 1:size(sigvecs,2)
        [~, maxpos] = (max(abs(sigvecs(peakreg:end,jj))));
        if sigvecs((peakreg + maxpos -1),jj) * peakstasign(ii) < 0
            sigvecs(:,jj) = -1* sigvecs(:,jj);
        end
    end
    % calculate nonlinearity for the eigen vectors
    sigvecsnorm = bsxfun(@times,sigvecs,1./sqrt(sum(sigvecs.^2,1)));
    siggen = res.runnoiseEn' * sigvecsnorm;
    stcout{ii}.sigvecs = sigvecs;
    stcout{ii}.sigvecsnorm = sigvecsnorm;
    for jj = 1:size(sigvecs,2)
        [nly, nlx] = rfroutines.getNonlinearity(siggen(:,jj),res.runnoisespkbins(ii,:), nlbin, toHzratio);        
        stcout{ii}.nly(jj,:) = nly;
        stcout{ii}.nlx(jj,:) = nlx;
         if sum(nly>0)>0
        stcout{ii}.predictions(jj,:) = rfroutines.getPredictionFromBinnedNonlinearity(res.frozenGenerators(ii,:),nlx,nly);
        [stcout{ii}.Rsq(jj),stcout{ii}.pearsonCoeff(jj),stcout{ii}.explVar(jj)] = ...
            calcRsqPearsonCoeff(res.frozenRates(ii,:),stcout{ii}.predictions(jj,:));
        cellrates = squeeze(res.trialrates(ii,:,:));       
        stcout{ii}.CCnorm(jj)    = rfroutines.calc_CCnorm(cellrates',stcout{ii}.predictions(jj,:)');
         end
    end
end
fprintf('Done! Took %2.2f s...\n', toc);
stcout = cell2mat(stcout);
end

%--------------------------------------------------------------------------------------------------%

%--------------------------------------------------------------------------------------------------%

function [corrstim1,corrstim2,combstim] = apply1stlevelNL(sta1,sta2, stim1,stim2, NLtype, para) %#ok

%if iscolumn(stim1), stim1 = transpose(stim1); end
filtData1 = filter2(sta1,stim1,'full') / para.resolution;
convstim1 = filtData1(1:end-(para.numSTA-1));

%if iscolumn(stim2), stim2 = transpose(stim2); end
filtData2 = filter2(sta2,stim2,'full') / para.resolution;
convstim2 = filtData2(1:end-(para.numSTA-1));

corrstim1 = applyNonlinearFilter(convstim1, NLtype);
corrstim2 = applyNonlinearFilter(convstim2, NLtype);

% summing each corrected nonlinearty
combstim = corrstim1 + corrstim2;

end

%--------------------------------------------------------------------------------------------------%

function [poisspktimes, poispsth] = poissonSpikeGenerator(predsignal, para,varargin)        %#ok
%
timeRes = 0.05;     % time bins in ms
repscaleRatio = round((1/(60/para.nblinks)*1000) / timeRes);

if iscolumn(predsignal), predsignal = transpose(predsignal); end

highrespsth = kron(predsignal,ones(1,repscaleRatio * para.nRepeats));       % replicating the predsignal
poisvals = poissrnd(highrespsth/((1/timeRes)*1e3));  % because the original predsignal is Hz and now is converted to ms
poisspk = reshape(poisvals,para.nRepeats,[]);

% get the spike timings
poisspktrials = cell(1,para.nRepeats);
for j = 1:para.nRepeats
    poisspktrials{j} = find(poisspk(j,:) == 1);
end
poisspktimes = CelltoMatUE(poisspktrials) * ((1/60 * para.nblinks)/repscaleRatio);  % in seconds
% generating psth from the poissonian rasters
poispsth = mean(reshape(mean(poisspk),repscaleRatio,length(predsignal)))*((1/timeRes)*1e3);

end

%--------------------------------------------------------------------------------------------------%

function plotFrozenNoiseMonochromatic(fndataall, para, clus, savingpath)

predtime = fndataall(1).frozentvec;
tvec = -(fndataall(1).statvec);%fliplr(linspace(0,para.numSTA*(1e3/(60/para.nblinks)),para.numSTA));
nstf = @(x)(num2str(round(x,2)));
msg =[];
for ii = 1: size(clus,1)
    fndata = fndataall(ii);
    nlmax = ceil(max(fndata.nly)/2)*2;
    if isnan(nlmax) ||(nlmax <= 0), nlmax = 10; end
    psmax = ceil(max([fndata.frozenpsth,fndata.pred])/2)*2;
    if isnan(psmax) ||(psmax <= 0), psmax = 10; end
    
    h = figure('position',[100 10 1700 850],'color',[1 1 1],'visible','off');
    
    subplot_tight(2,4,1,0.05)
    plot(tvec,fndata.stanorm,'color',rgb('red'),'linewidth',2);
    axis([0 0.5 -0.75 0.75]); set(gca,'xtick',0:0.250:1,'ytick',-2:0.5:2);   axis square;
    ylabel('filter strength');     title('STA');      %   xlabel('time (sec)');
    
    subplot_tight(2,4,5,0.05);
    plot(fndata.nlx,fndata.nly,'-o','color',rgb('red'),'markerfacecolor',rgb('salmon'),'linewidth',2);
    axis([-3 3 0 nlmax]);      set(gca,'xtick',-6:1.5:6,'ytick',0:nlmax/2:nlmax);     axis square;
    ylabel('Hz');   xlabel('input');       title('Nonlinearity');
    
    subplot_tight(2,4,[2 4],0.05)
    rasterPlotter(fndata.frozenspikes,0.4,rgb('royalblue'),0.25);
    axis([min(predtime) max(predtime) 0 ceil(size(fndata.frozenspikes,2)/2)*2]);
    ax = gca;   ax.YTick = 0:(ceil(size(fndata.frozenspikes,2)/2)*2)/2:ceil(size(fndata.frozenspikes,2)/2)*2;
    ax.XColor = 'w';
    ax.TickLength = [0.0025 0.01];     ylabel('trials');     title('raster plot from frozen noise');
    
    subplot_tight(2,4,[6 8],0.05)
    plot(predtime, fndata.frozenpsth,'color',rgb('royalblue'),'linewidth',1);
    hold on
    plot(predtime, fndata.pred,'color',rgb('red'),'linewidth',1.5);
    axis([min(predtime) max(predtime) 0 psmax]);
    box off;    ax = gca;    ax.XTick = 0:max(round(predtime,1))/2:max(round(predtime,1));
    xlabel('time (sec)');    ax.TickLength = [0.0025 0.01];
    ylabel('Hz');  title(['prediction vs response, Rsq: ',nstf(fndata.predstat.pearCorrSqrt)]);
    legend('psth','prediction','numcolumns',2);     legend('boxoff');
    % saving data and plot
    filename = generateRGCname('Frozen Noise',clus(ii,:),savingpath);
    suptitle_mod(h,filename,3);
    
    savepngFast(h,savingpath,filename);
    close(h);
    fndata.para = para;
    save([savingpath,'/fna_data/',filename,'.mat'],'-struct','fndata');
    %disp(['Analysis for Cell ',num2str(clus(ii,1)),', Cluster ',num2str(clus(ii,2)),' is... wait for it...done!!!']);
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', clus(ii,1),clus(ii,2));
    fprintf(msg);
end

end

%--------------------------------------------------------------------------------------------------%

function plotFrozenNoiseDichromatic(fndataall, ciff, para, clus, savingpath)

predtime = fndataall(1).frozentvec;
%linspace(para.numSTA *para.nblinks/60 ,para.FrozenFrames* para.nblinks/60, para.FrozenFrames-para.numSTA+1);
tvec = -(fndataall(1).statvec);%fliplr(linspace(0,para.numSTA*(1e3/(60/para.nblinks)),para.numSTA));
nstf = @(x)(num2str(round(x,2)));

for ii = 1: size(clus,1)
    fndata = fndataall(ii);
    nlmax = ceil(max([fndata.nlygreen,fndata.nlyuv,fndata.nlyjoined])/2)*2;
    if isnan(nlmax) ||(nlmax <= 0), nlmax = 10; end
    psmax = ceil(max([fndata.frozenpsth,fndata.predjoined,fndata.predgreen,fndata.preduv])/2)*2;
    if isnan(psmax) ||(psmax <= 0), psmax = 10; end
    
    h = figure('position',[100 10 1700 1050],'color',[1 1 1],'visible','off');
    
    subplot_tight(4,5,1,0.05)
    plot(tvec,fndata.stagreennorm,'color',rgb('green 3'),'linewidth',1);
    axis([0 0.5 -0.75 0.75]); set(gca,'xtick',0:0.250:1,'ytick',-2:0.5:2);   axis square;
    ylabel('filter strength');  title('STA green');
    
    subplot_tight(4,5,2,0.05);
    plot(fndata.nlxgreen,fndata.nlygreen,'color',rgb('green 3'),'linewidth',1);
    axis([-3 3 0 nlmax]);      set(gca,'xtick',-6:1.5:6,'ytick',0:nlmax/2:nlmax);     axis square;
    ax = gca;   ax.Position(1) = ax.Position(1)-0.025;
    ylabel('Hz');  title('Nonlinearity green');
    
    subplot_tight(4,5,[3 5],[0.05,0.01])
    plot(predtime, fndata.frozenpsth,'color',rgb('dimgray'));
    hold on
    plot(predtime, fndata.predgreen,'color',rgb('green 3'));
    axis([min(predtime) max(predtime) 0 psmax]);
    box off;    ax = gca;   ax.XColor = 'w';    ax.XLabel = []; ax.YTick = 0:psmax/2:psmax;
    ax.TickLength = [0.0025 0.01];   ax.Position(1) = ax.Position(1)-0.01;
    ylabel('Hz');  title(['prediction green vs response, pCoeff: ',nstf(fndata.predgreenstat.pearCorrSqrt)]);
    
    subplot_tight(4,5,6,0.05)
    plot(tvec,fndata.stauvnorm,'color',rgb('blueviolet'),'linewidth',1);
    axis([0 0.5 -0.75 0.75]); set(gca,'xtick',0:0.25:1,'ytick',-2:0.5:2);     axis square;
    ylabel('filter strength');  title('STA UV');
    
    subplot_tight(4,5,7,0.05);
    plot(fndata.nlxuv,fndata.nlyuv,'color',rgb('blueviolet'),'linewidth',1);
    axis([-3 3 0 nlmax]);      set(gca,'xtick',-6:1.5:6,'ytick',0:nlmax/2:nlmax);   axis square;
    ax = gca;   ax.Position(1) = ax.Position(1)-0.025;
    ylabel('Hz');  title('Nonlinearity UV');
    
    subplot_tight(4,5,[8 10],[0.05,0.01])
    plot(predtime, fndata.frozenpsth,'color',rgb('dimgray'));
    hold on
    plot(predtime, fndata.preduv,'color',rgb('blueviolet'));
    axis([min(predtime) max(predtime) 0 psmax]);
    box off;    ax = gca;   ax.XColor = 'w';    ax.XLabel = []; ax.YTick = 0:psmax/2:psmax;
    ax.TickLength = [0.0025 0.01];           ax.Position(1) = ax.Position(1)-0.01;
    ylabel('Hz');  title(['prediction UV vs response, pCoeff: ',nstf(fndata.preduvstat.pearCorrSqrt)]);
    
    subplot_tight(4,5,11,0.05)
    plot(tvec,fndata.stabothnormgr,'color',rgb('green 3'),'linewidth',1);
    hold on
    plot(tvec,fndata.stabothnormuv,'color',rgb('blueviolet'),'linewidth',1);
    axis([0 0.5 -0.75 0.75]); set(gca,'xtick',0:0.25:1,'ytick',-2:0.5:2);        axis square;
    ylabel('filter strength'); title('STA green & UV');
    
    subplot_tight(4,5,12,0.05);
    plot(fndata.nlxjoined,fndata.nlyjoined,'color',rgb('red'),'linewidth',1);
    axis([-3 3 0 nlmax]);      set(gca,'xtick',-6:1.5:6,'ytick',0:nlmax/2:nlmax);    axis square;
    ax = gca;   ax.Position(1) = ax.Position(1)-0.025;
    ylabel('Hz');  xlabel('input');  title('Nonlinearity joined (green + UV)');
    
    subplot_tight(4,5,[13 15],[0.05,0.001])
    plot(predtime, fndata.frozenpsth,'color',rgb('dimgray'));
    hold on
    plot(predtime, fndata.predjoined,'color',rgb('red'));
    axis([min(predtime) max(predtime) 0 psmax]);
    box off;    ax = gca;   ax.XColor = 'w';    ax.XLabel = []; ax.YTick = 0:psmax/2:psmax;
    ax.TickLength = [0.0025 0.01];      ax.Position(1) = ax.Position(1)-0.01;
    ylabel('Hz');  title(['prediction joined vs response, pCoeff: ',nstf(fndata.predjoinedstat.pearCorrSqrt)]);
    
    subplot_tight(4,5,16,0.05)
    [~,~,cip] = plotColorIntegration(ciff(ii).stimpsth, ciff(ii).pfrpsth,'markersize',4,'linewidth',1,'showylabel',0,'facealpha',0.2);
    xticklabels(1:1:11);    ylabel('spike differences');  xlabel('contrast combinations'); title('chromatic integration curves');
    
    subplot_tight(4,5,[18 20],[0.05,0.01])
    rasterPlotter(fndata.frozenspikes,0.4,rgb('dimgray'),0.25);
    axis([min(predtime) max(predtime) 0 ceil(size(fndata.frozenspikes,2)/2)*2]);
    ax = gca;   ax.YTick = 0:(ceil(size(fndata.frozenspikes,2)/2)*2)/2:ceil(size(fndata.frozenspikes,2)/2)*2;
    ax.TickLength = [0.0025 0.01];      ax.Position(1) = ax.Position(1)-0.01;   ax.XTick = 0:max(round(predtime,1))/2:max(round(predtime,1));
    xlabel('time (sec)'); ylabel('trials');     title('raster plot from frozen noise');
    
    subplot_tight(4,5,17,0.05)
    text(-.08,0.9,{['Pearson Coeff joined STAs: ',nstf(fndata.predjoinedstat.pearCorrSqrt)]; ['Nonlinearity index: ',nstf(ciff(ii).indices.nlindex)];...
        'Nonlinearity type :';['     ' ciff(ii).ciclasstype];['Pearson Coeff green STA: ',nstf(fndata.predgreenstat.pearCorrSqrt)];...
        ['Pearson Coeff UV STA: ',nstf(fndata.preduvstat.pearCorrSqrt)];['Rsq joined STAs: ',nstf(fndata.predjoinedstat.Rsqr)]; ...
        ['Rsq green STA: ',nstf(fndata.predgreenstat.Rsqr)]; ['Rsq UV STA: ',nstf(fndata.preduvstat.Rsqr)]},'horiz','left','vert','top');
    axis off;
    l = legend([cip.plot1,cip.plot2],'Green-ON, UV-OFF','Green-OFF, UV-ON');
    l.Position = [l.Position(1)+0.14,l.Position(2)-0.17, l.Position(3), l.Position(4)];     l.Box = 'off';
    % saving data and plot
    filename = generateRGCname('Frozen Noise',clus(ii,:),savingpath);
    suptitle_mod(h,filename,3);
    
    savepngFast(h,savingpath,filename);
    close(h);
    fndata.para = para;
    save([savingpath,'/fna_data/',filename,'.mat'],'-struct','fndata');
    disp(['Analysis for Cell ',num2str(clus(ii,1)),', Cluster ',num2str(clus(ii,2)),' is... wait for it...done!!!']);
    
end

end

