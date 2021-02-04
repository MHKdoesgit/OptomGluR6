

function varargout = Analyze_Chirp(datapath, varargin)
%
%%% Analyze_Chirp %%%
%
%
% This function analyze chirp stimulus. The stimulus is created by baden et al,
% and it has on-off step, temporal frequency sweep and contrast sweep part.
% This function additionally measures the frequency senitivity and contrast
% sensitiviy of the ganglion cells.
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   chripdata : a structure containing rasters values, psth values and
%             output of spikes differences analysis.
%   Plot : This function will also plot the rasters and psth and the other
%          relevant plots.
%
% written by Mohammad, started ages ago but finished on 18.07.2020, after
% possibly 1-2 years of postponing.

totaltime =tic;
if (nargin < 1),  datapath = uigetdir();       end

[ftimes, spikes,para, clus, savingpath] = chirp_parameters(datapath,varargin);
% first analyze and save everything
res = analyze_chrip_stim(ftimes,spikes, clus, para, savingpath);
% then plot all togetger
plot_chirp_stimulus(res, clus, savingpath);

sound(struct2array(load('gong.mat','y')));
disp(seconds2human (toc(totaltime)));
varargout{1} = res;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [ft, spikes,para, clus, savingpath] = chirp_parameters(dp,varargin)

[thisExp, savingpath] = loadRawData(dp,'chirp','Chirp_Analysis');
% loading clusters and frametimings
clus = thisExp.clusters;
ft = thisExp.ftimes;

para.onoff = thisExp.stimPara.onoffstep;
para.onoff.numframes = 4;
para.freqsweep = thisExp.stimPara.freqsweep;
para.contrastsweep = thisExp.stimPara.ctrsweep;
%para.colorFlag = thisExp.stimPara.color;
para.fps = thisExp.info.screen.refreshrate;
if datenum(thisExp.date) < datenum('01-Jan-2021')
    para.stimfps = 60;
else
    para.stimfps = thisExp.info.screen.refreshrate;
end
para.binlength = 10./1e3;
if thisExp.stimPara.color, para.coneisolating = thisExp.stimPara.coneisolating; end
%para.Nrepeats = thisExp.stimPara.Nrepeats;

%para.chagetoEulerdate = '01-Nov-2017';
para.Nblinks = 2; % defiend as constant in the stimulus program
if datenum(thisExp.date) < datenum('01-Nov-2017')
    
    para.trialLength = para.onoff.numframes + para.freqsweep.duration/2 + para.contrastsweep.duration/2;

elseif datenum(thisExp.date) < datenum('01-Jan-2021') && datenum(thisExp.date) > datenum('01-Nov-2017')
    
    para.Nrepeats = length(find(diff(ft)> (para.onoff.preframes-2) /para.fps));
    stfr = (find(diff(ft)>((para.Nblinks+1)/para.fps)));
    onetriallen = ceil(length(stfr)/para.Nrepeats);
    if length(stfr) ~= onetriallen*para.Nrepeats
        stfr = (find(diff(ft)>((para.Nblinks+2)/para.fps)));
        %             difflen = onetriallen*para.Nrepeats - length(stfr);
        %             stfr = [stfr; nan(difflen,1)];
    end
    changepoints = reshape(stfr,[],para.Nrepeats)'; % anything bigger than 2blinks is special
    chpts  = [changepoints(:,1:4),changepoints(:,4)+1, changepoints(:,5), changepoints(:,5)+1 ,...
        [changepoints(2:end,1)-1; length(ft)]]; 
    % these are all the stimulus change points ===>
    %|1--gray--|2--on--|3--off--|4--gray--|5--frqsweeps--|245--gray--|246--contsteps--|515
    % now remove the third frame, between the on and off and find the frame
    % for the gray screen before the end!
    graybeforeend = chpts(:,end-1)+(1/(para.Nblinks/para.fps)*(para.contrastsweep.duration/para.fps));
    %graybeforeend = chpts(:,end)+1;
    %para.trials = [chpts(:,1:2),chpts(:,4:end-1),graybeforeend, chpts(:,end)];
    para.trials = [chpts(:,1:end-1),graybeforeend,chpts(:,end)];
    para.ftchanges = chpts;
else
    %|1--gray--|2--on--|3--off--|4--gray--|5--frqsweeps--|6--gray--|7--contsteps--|8--gray--|9--end--
    numpulsepertrial = 9;
    para.Nrepeats = floor(length(ft)/numpulsepertrial);
    stfr = ft(1:numpulsepertrial*para.Nrepeats);
    changepoints = reshape(1:numel(stfr),[],para.Nrepeats)';
    chpts = changepoints-transpose(0:para.Nrepeats-1);    
    para.trials = chpts;
    para.ftchanges = chpts;
end

tvec = (ft(para.trials(1,:))-ft(para.trials(1,1)));
%tvec = [(tvec-2),tvec(end)];
%tvec(tvec<0) = 0;
%tvec = [tvec(1:4),tvec(4)+2,tvec(5),tvec(5)+2,tvec(6)-2,tvec(6:end)];
para.changepoints = tvec';

% for frequency sweeps
% from the stimulus program
% phase = t * (loFreq + ((hiFreq - loFreq) / incDuration) * t / 2);
nstimbins = para.freqsweep.duration/para.fps / para.binlength;
para.freqsweep.tvec = linspace(0, para.freqsweep.duration/para.stimfps, nstimbins);
para.freqsweep.freq = linspace(para.freqsweep.low_freq,para.freqsweep.high_freq, nstimbins);
para.freqsweep.stim = sin(2*pi.* para.freqsweep.freq.* (para.freqsweep.tvec/2));
[para.freqsweep.peaks,para.freqsweep.peaklocs] = findpeaks(para.freqsweep.stim);
para.freqsweep.fitopts = optimset('Display','off','TolX', 1e3*eps,'TolFun', 1e3*eps,'MaxFunEvals',2e4,'MaxIter',2e4);
para.freqsweep.predFun=@(p,x) (p(1)*(1-exp(-(x.^2)./(2*(p(3).^2))))-(p(2)*(1-exp(-(x.^2)./(2*(p(4).^2))))));

% for contrast sweeps
% from stimulus program
% contrast = ctrMaxContrast * ((double)t / ctrDuration) * sin(2 * Math::PI * freq*t / fps);
para.contrastsweep.maxcontrast = para.contrastsweep.contrast;
nstimbins = para.contrastsweep.duration/para.fps / para.binlength;
para.contrastsweep.tvec = linspace(0, para.contrastsweep.duration/para.stimfps,nstimbins);
freq = repmat(para.contrastsweep.freq, 1, nstimbins);
para.contrastsweep.stim = para.contrastsweep.maxcontrast * sin(2*pi.*freq.*para.contrastsweep.tvec )...
    .*(para.contrastsweep.tvec /(para.contrastsweep.duration/para.stimfps));
[para.contrastsweep.peakspos,para.contrastsweep.peaklocspos]= findpeaks(para.contrastsweep.stim ); % positive peaks
[yneg,para.contrastsweep.peaklocsneg]= findpeaks(-para.contrastsweep.stim ); % negative peaks
para.contrastsweep.peaksneg = -yneg;
para.contrastsweep.contrast = linspace(0,para.contrastsweep.maxcontrast,nstimbins);
para.contrastsweep.fitopts = optimset('Display','off','TolX', 1e3*eps,'TolFun', 1e3*eps,'MaxFunEvals',...
    2e4,'MaxIter',2e4,'Jacobian', 'on');
para.samplingrate = thisExp.info.samplingrate;
if isfield(thisExp,'sortinginfo')
    para.sortinfo = thisExp.sortinginfo;
else
    para.sortinfo = [];
end
para.date = thisExp.date;
para.expnumber = thisExp.stimPara.expnumber;

spikes = thisExp.spiketimes;

end

%--------------------------------------------------------------------------------------------------%

function res = analyze_chrip_stim(ft, spikes, clus, para, savingpath, varargin)

tic;
stimft = ft(para.trials);
nbins = (round(stimft(1,1:end)-stimft(1,1))/para.binlength);
nbins = [diff(nbins),sum(diff(nbins))];
nbintr = nbins(end)+1;      nbins = nbins(1:end-1)+1;
binsdur = diff(para.changepoints);
%nbins = [nbins(1), sum(nbins(2:3)), nbins(4:end)];

para.chirpstim = [zeros(1,nbins(1)-1),ones(1,(nbins(2)-1)),-1*ones(1,(nbins(3)-1)),zeros(1,nbins(3)-1),...
    para.freqsweep.stim,zeros(1,nbins(6)-1),para.contrastsweep.stim,zeros(1,nbins(8)-1)];
para.nbins = [nbins, nbintr];

ptimes.grfirst = linspace(para.changepoints(1),para.changepoints(2),nbins(1)-1);
ptimes.onpsth = linspace(para.changepoints(2),para.changepoints(3),nbins(2)-1);
ptimes.offpsth = linspace(para.changepoints(3),para.changepoints(4),nbins(3)-1);
ptimes.grfrq = linspace(para.changepoints(4),para.changepoints(5),nbins(4)-1);
ptimes.frqpsth = linspace(para.changepoints(5),para.changepoints(6),nbins(5)-1);
ptimes.grcont = linspace(para.changepoints(6),para.changepoints(7),nbins(6)-1);
ptimes.contpsth = linspace(para.changepoints(7),para.changepoints(8),nbins(7)-1);
ptimes.grend = linspace(para.changepoints(8),para.changepoints(9),nbins(8)-1);
ptimes.psth = linspace(para.changepoints(1),para.changepoints(end),nbintr-1);

prt =  zeros(size(clus,1),nbintr, para.Nrepeats);
grfirst = zeros(size(clus,1), nbins(1)-1);
onpsth = zeros(size(clus,1), nbins(2)-1);
offpsth = zeros(size(clus,1), nbins(3)-1);
grfrq = zeros(size(clus,1), nbins(4)-1);
frqpsth = zeros(size(clus,1), nbins(5)-1);
grcont = zeros(size(clus,1), nbins(6)-1);
contpsth = zeros(size(clus,1), nbins(7)-1);
grend = zeros(size(clus,1), nbins(8)-1);
rasters = cell(size(clus,1), size(stimft,2));
Rsqs = nan(size(clus,1), size(stimft,2));
% [frqpeaks, frqplocs, frqfit] = deal(zeros(size(clus,1), length(para.freqsweep.peaks)));
% [contpeaks, contplocs, contfit] = deal(zeros(size(clus,1), length(para.contrastsweep.peakspos)));
% [frqlagresp, contlagresp] = deal(zeros(size(clus,1),1));
% [frqfitparams, contfitparams] = deal(zeros(size(clus,1),4));

for ii = 1:size(clus,1)
    
    spk = spikes{ii};
    
    rastrial = cell(para.Nrepeats,1);
    rasparts = cell(para.Nrepeats,size(para.trials,2)-1);
    psthparts = cell(size(para.trials,2)-1,1);
    for jj = 1:para.Nrepeats
        rastrial{jj} = spk(and(spk > stimft(jj,1), spk <= stimft(jj,end))) - stimft(jj,1);
        if isempty(rastrial{jj}), continue;end
        prt(ii,:,jj) = histc(rastrial{jj}',linspace(para.changepoints(1),para.changepoints(end),nbintr)); %#ok
        for kk = 1:size(stimft,2)-1
            rasparts{jj,kk} = spk(and(spk > stimft(jj,kk), spk <= stimft(jj,kk+1))) - stimft(jj,kk);
            
        end
    end
    rasters{ii,end} = CelltoMatUE(rastrial);
    
    for kk = 1:size(stimft,2)-1
        rasters{ii,kk} = CelltoMatUE(rasparts(:,kk));
        if size(rasters{ii,kk},2)==1
            psthparts{kk} = mean(histc(rasters{ii,kk}',linspace(0,binsdur(kk),nbins(kk))),1)'; %#ok, has to be histc for 2d conversion
        else
            psthparts{kk} = mean(histc(rasters{ii,kk}',linspace(0,binsdur(kk),nbins(kk))),2); %#ok
        end
        psthparts{kk} = psthparts{kk}(1:end-1)*(1/para.binlength);
        %rsqras = rasters{ii,kk};  if size(rsqras,2)==1, rsqras = [rsqras, nan(para.Nrepeats,1)]; end
        rsqras = rasters{ii,kk};        rsqras(isnan(rsqras)) = 0;
        Rsqs(ii,kk) = rfroutines.imageTrialRsq(reshape(rsqras,1,para.Nrepeats,[]));
    end
    rsqras = rasters{ii,end};        rsqras(isnan(rsqras)) = 0;
    Rsqs(ii,end) = rfroutines.imageTrialRsq(reshape(rsqras,1,para.Nrepeats,[]));
    
    grfirst(ii,:) = psthparts{1};    onpsth(ii,:) = psthparts{2};   offpsth(ii,:) = psthparts{3};
    grfrq(ii,:) = psthparts{4};      frqpsth(ii,:) = psthparts{5};
    grcont(ii,:) = psthparts{6};     contpsth(ii,:) = psthparts{7};
    grend(ii,:) = psthparts{8};
    
    %     if ~isempty(psthparts{5})
    %         [frqpeaks(ii,:), frqplocs(ii,:), frqlagresp(ii)] = chirppeaks(psthparts{5}, para.freqsweep.stim);
    %
    %         [frqfit(ii,:), frqfitparams(ii,:)] = fitDoG(para.freqsweep.freq(frqplocs(ii,:)),frqpeaks(ii,:), para);
    %     end
    %     if ~isempty(psthparts{7})
    %         [contpeaks(ii,:), contplocs(ii,:), contlagresp(ii)] = chirppeaks(psthparts{7}, para.contrastsweep.stim);
    %
    %         [contfit(ii,:), contfitparams(ii,:)] = fitsigm(para.contrastsweep.contrast(contplocs(ii,:)),contpeaks(ii,:), para);
    %     end
end

prt = prt(:,1:nbintr-1,:);
ratesOdd=mean(prt(:,:,1:2:end),3);
ratesEven=mean(prt(:,:,2:2:end),3);
qualityRsq=1-sum((ratesOdd-ratesEven).^2,2)./sum((ratesOdd-repmat(mean(ratesOdd,2),[1 size(ratesOdd,2)])).^2,2);
toc
psth = mean(prt, 3) * (1/para.binlength);

res.rasters = rasters;
res.psth = psth;
res.qualityRsq = qualityRsq;
res.rsqall = Rsqs;
res.onoff.onpsth = onpsth;
res.onoff.offpsth = offpsth;
% freq sweeps
res.freqsweep.psth = frqpsth;
res.freqsweep.rasterindex = 5;
%frqpeakfit = frqsweepAmp(frqpsth, para);
res.freqsweep = struct2struct(res.freqsweep, frqsweepAmp(frqpsth, para));

% res.freqsweep.peaks = frqpeaks;
% res.freqsweep.peaklocs = frqplocs;
% res.freqsweep.timetoresp = frqlagresp;
% res.freqsweep.fit = frqfit;
% res.freqsweep.fitparams = frqfitparams;
% contrast sweeps
res.contsweep.psth = contpsth;
res.contsweep.rasterindex = 7;
res.contsweep = struct2struct(res.contsweep, contsweepAmp(contpsth,  para));
% res.contsweep.peaks = contpeaks;
% res.contsweep.peaklocs = contplocs;
% res.contsweep.timetoresp = contlagresp;
% res.contsweep.fit = contfit;
% res.contsweep.fitparams = contfitparams;
% preframes
res.preframe.onoff = grfirst;
res.preframe.freqsweep = grfrq;
res.preframe.contsweep = grcont;
res.preframe.final = grend;
res.psthtimes = ptimes;
res.clusters = clus;
res.para = para;

chripdata = res;
filename = [num2str(para.expnumber,'%02d'),'-Chirpstimulus_analysis_for_experiment_on_',para.date,'.mat'];
save([savingpath,filename],'-struct','chripdata');


%--------------------------------------------------------------------------------------------------
% bd=load('baden_data_repurposed.mat');
% kernelTau=0.5664;
% kernelConst=0.16;
% kernelTmax=4;
% tKernel=0:para.Nblinks/para.fps:kernelTmax;
% caKernel=(1-kernelConst)*exp(-tKernel/kernelTau)+kernelConst;
% caTimes=bd.chirp_time;
%
% caSignals = getCalciumTracesNew(prt(:,1:2300,:), ptimes(1:2300), caKernel,caTimes);
% cMat=corr(caSignals',bd.chirp_group_mean', 'type','Pearson');
% res.caSignals=caSignals;
% res.caTimes=caTimes;
% res.corrMat=cMat;

% calciumevent = 0.15+[zeros(1,2000),(exp((0:-0.001:(-5+0.001))/0.57))];
% casig = zeros(para.Nrepeats, para.changepoints(end)*1e3);
% for jj = 1:para.Nrepeats
%
%     psthpertrial = histc(rasters.trialraster(jj,:),0:0.001:para.changepoints(end));
%     convsig = conv(psthpertrial, calciumevent);
%     casig(jj,:) = convsig(length(calciumevent)/2:end-(length(calciumevent)/2+1));
% end
%
% ca2signal.trialsignals = casig;
% ca2signal.meanca2signal = mean(casig,1);
% ca2signal.ca2event = calciumevent;

% p = mean(histc(rasters.trialraster',0:0.001:32),2);
% cs = conv(p,calciumevent);

end

%--------------------------------------------------------------------------------------------------%

function out = frqsweepAmp(frqpsth, para)

disp('analyzing frequency sweeps...');
[~,fl] =findpeaks(-abs(para.freqsweep.stim));
zerolocs = [1,fl];
periodedges = zerolocs(1:2:end);
periodedgesend = [periodedges(2:end), length(para.freqsweep.stim)];

stimidx = arrayfun(@(x,y)x:y,periodedges, periodedgesend,'un',0);

frqloc = para.freqsweep.freq(zerolocs(2:2:end));

maxeachfreq = nan(size(frqpsth,1), length(periodedges));
%respamps = frqpsth ./ max(frqpsth,[],2);

for jj = 1:length(periodedges)
    maxeachfreq(:,jj) = max(frqpsth(:,stimidx{jj}),[],2);
end

xvals = linspace(para.freqsweep.low_freq,para.freqsweep.high_freq,100);
uparams = nan(size(frqpsth,1), 4);
mdlresp = zeros(size(frqpsth,1),100);
for ii = 1:size(frqpsth,1)
    uparams(ii,:) = fitGenLogisticToSpikes(frqloc, maxeachfreq(ii,:));
    mdlresp(ii,:) = logistic4(uparams(ii,:), xvals);
end

out.maxeachfreq = maxeachfreq;
out.freqs = frqloc;
out.fitparams = uparams;
out.modelresp = mdlresp;
out.fitxax = xvals;
out.periodedges = [periodedges;periodedgesend]';

end

%--------------------------------------------------------------------------------------------------%

function out = contsweepAmp(contpsth,  para)

disp('analyzing contrast sweeps...')
[~,fl] =findpeaks(-abs(para.contrastsweep.stim));
zerolocs = [1,fl];
periodedges = zerolocs(1:2:end);
periodedgesend = [periodedges(2:end), length(para.contrastsweep.stim)];
stimidx = arrayfun(@(x,y)x:y,periodedges,periodedgesend,'un',0);

contloc = para.contrastsweep.contrast(zerolocs(2:2:end));

maxeachperiod = nan(size(contpsth,1), length(periodedges));
%respamps = contpsth ./ max(contpsth,[],2);

for jj = 1:length(periodedges)
    maxeachperiod(:,jj) = max(contpsth(:,stimidx{jj}),[],2);
end

xvals = linspace(min(para.contrastsweep.contrast),para.contrastsweep.maxcontrast,100);
uparams = nan(size(contpsth,1), 4);
mdlresp = zeros(size(contpsth,1),100);
for ii = 1:size(contpsth,1)
    uparams(ii,:) = fitGenLogisticToSpikes(contloc, maxeachperiod(ii,:));
    mdlresp(ii,:) = logistic4(uparams(ii,:), xvals);
end

out.maxeachperiod = maxeachperiod;
out.contrasts = contloc;
out.fitparams = uparams;
out.modelresp = mdlresp;
out.fitxax = xvals;
out.periodedges = [periodedges;periodedgesend]';

% plot(contloc,maxeachperiod(ii,:),'o')
% hold on
% plot(xvals,mdlresp ./ max(mdlresp))
%
end

%--------------------------------------------------------------------------------------------------%

function [yp, xp, lagresp, xc, lag] = chirppeaks(psth, stim) %#ok

[xc,lag] = xcorr(psth,stim);
[~,lagresp] = max(xc);
lagresp = lag(lagresp);

[~, nxp] = findpeaks(-stim,lagresp+(1:length(stim)));

if sum(nxp>length(psth)) > 1 || lagresp < 0
    [~, nxp] = findpeaks(-stim,5+(1:length(stim)));
end
if nxp(end)> length(psth), nxp(end) = length(psth)+1; end
nxp = [1,nxp];
idx = @(x1,x2)(nxp(x1)+1: nxp(x2)-1);
[yp,xp] = deal(zeros(size(nxp,2)-1,1));
for kk = 1:length(nxp)-1
    kkidx = idx(kk,kk+1);
    if isempty(psth), continue; end
    if isempty(kkidx) && kk==length(nxp)-1, kkidx = length(psth)-3:length(psth); end
    while length(kkidx)<3, kkidx= kkidx(1)-1:kkidx(end); end
    [ypp, xpp] = findpeaks(psth(kkidx),'Npeak',1,'SortStr','descend') ;
    if isempty(xpp)
        xpp = floor(median(kkidx))-nxp(kk);
        ypp = max(psth(xpp+nxp(kk)));
    end
    yp(kk) = ypp;
    xp(kk) = kkidx(1) + xpp;
end

end

%--------------------------------------------------------------------------------------------------%

function [ypred, params] = fitsigm(xx,yy, para) %#ok

guess = [yy(1) max(yy) 15 1];
try
    params = lsqcurvefit(@(p,x) logistic4(p,x), guess, xx,yy', [0 0 0 0], [5 1000 300 Inf], ...
        para.contrastsweep.fitopts);
catch ME
    disp(ME.message);
    params = lsqcurvefit(@(p,x) logistic4(p,x), guess, xx,yy');
end
ypred = logistic4(params,xx);

end

%--------------------------------------------------------------------------------------------------%

function [ypred, params] = fitDoG(xx,yy, para) %#ok

guess = [max(yy) max(yy)*0.5 1 1]; % 250 for center and 750 for surround
% predFun=@(p,x) (p(1)*(1-exp(-(x.^2)./(2*(p(3).^2))))-(p(2)*(1-exp(-(x.^2)./(2*(p(4).^2))))));
try
    params = lsqcurvefit(para.freqsweep.predFun, guess, xx,yy, [0 0 0 0], ...
        [max(yy)*2 max(yy)*2 100 Inf], para.freqsweep.fitopts);
    
catch ME
    disp(ME.message);
    params = lsqcurvefit(para.freqsweep.predFun, guess, xx,yy);
end

ypred = para.freqsweep.predFun(params, xx);

end

%--------------------------------------------------------------------------------------------------%

function [x,y] = raslin(r)

spkheight = 0.45;
[x,y] = deal(cell(size(r)));
for ii = 1: size(r,1)
    for jj = 1: size(r,2)
        rasresp = r{ii,jj};     %CelltoMatUE(r(jj,:));
        spkloc = repmat(1:size(rasresp,1),1,size(rasresp,2));
        spkoneline = [rasresp(:),rasresp(:),nan(numel(rasresp(:)),1)]';
        tboneline = [spkloc(:)-spkheight,spkloc(:)+spkheight, nan(numel(spkloc(:)),1)]' ; % + rasloc(jj)
        
        x{ii,jj} = spkoneline(:);
        y{ii,jj} = tboneline(:);
    end
end

end

%--------------------------------------------------------------------------------------------------%

function plot_chirp_stimulus(res, clus, savingpath)

[x,y] = raslin(res.rasters);
p = res.para;
tt = 0.028;
offcol = rgb('royalblue');
oncol = rgb('crimson');
stimcol = rgb('red');
msg = [];

for ii = 1:size(clus,1)
    
    h = figure('pos',[200 10 1600 1000],'color','w','vis','off');
    % plot all rasters
    subplot_tight(4,1,1,tt)
    line(x{ii,end},y{ii,end},'color','k','linewidth',0.25);
    hold on;
    imagesc(linspace(0,p.changepoints(end),length(p.chirpstim)),p.Nrepeats+2:p.Nrepeats+2,repmat(p.chirpstim,2,1));
    colormap(gray);
    plot(linspace(0,p.changepoints(end),length(p.chirpstim)),p.chirpstim/2 + p.Nrepeats+4.3,'color',stimcol);
    axis([0 round(p.changepoints(end)), 0 p.Nrepeats+5]);
    ax = gca;       ax.XColor = 'none';     yticks(0:(p.Nrepeats+1)/2:p.Nrepeats+1);
    ax.TickLength= [0.002 0.002];           pbaspect([8 1 1]);      ylabel('trials');
    
    % plot psth
    subplot_tight(4,1,2,tt)
    line(res.psthtimes.psth,res.psth(ii,:),'color','k');
    xlim([0 round(p.changepoints(end))]);
    xticks(round(p.changepoints));      box off;
    ax = gca;       ax.TickLength= [0.002 0.002];
    line([p.changepoints(2:end);p.changepoints(2:end)],[0;ax.YLim(2)],'color','r','linestyle','--');
    pbaspect([8 1 1]);              xlabel('time (s)');         ylabel('rate');
    title(['Rsq overall: ',num2str(round(res.rsqall(ii,end),2))],'FontSize',9)
    %hold on;
    %plot(res.psthtimes.psth,smoothdata(res.psth(ii,:)),'r');
    
    % plot on-off steps rasters
    subplot_tight(4,4,9,tt)
    line(x{ii,1}, y{ii,1},'color',0.5* [1 1 1]);
    line(p.changepoints(2)+ x{ii,2}, y{ii,2},'color',oncol);
    line(p.changepoints(3)+ x{ii,3}, y{ii,3},'color',offcol);
    line(p.changepoints(4)+ x{ii,4}, y{ii,4},'color',0.5* [1 1 1]);
    xline(round(p.changepoints(2)),'--','on','color',oncol,'LabelHorizontalAlignment','left');
    xline(round(p.changepoints(3)),'--','off','color',offcol,'LabelHorizontalAlignment','left');
    xline(round(p.changepoints(4)),'--','preframe','color',0.5* [1 1 1],'LabelHorizontalAlignment','left');
    axis([0 round(p.changepoints(5)) 0 p.Nrepeats+1]);
    xticks(round(p.changepoints(1:5)));
    yticks(0:(p.Nrepeats+1)/2:p.Nrepeats+1);
    pbaspect([2 1 1]);           ax = gca;       ax.XColor = 'none';        ylabel('trials');
    title(['Rsq on: ',num2str(round(res.rsqall(ii,2),2)),', off: ',num2str(round(res.rsqall(ii,3),2))],'FontSize',9)
    
    % rasters from frequency sweeps
    subplot_tight(4,4,10,tt)
    line(x{ii,5}, y{ii,5},'color','k');
    line(res.psthtimes.frqpsth-res.psthtimes.frqpsth(1) , p.freqsweep.stim + p.Nrepeats+2,'color',stimcol);
    xticks(0:p.freqsweep.duration/p.fps/4:p.freqsweep.duration/p.fps);
    yticks(0:(p.Nrepeats+1)/2:p.Nrepeats+1);
    pbaspect([2 1 1]);           ax = gca;       ax.XColor = 'none';        ylabel('trials');
    title(['Rsq frequency sweep: ',num2str(round(res.rsqall(ii,5),2))],'FontSize',9)
    
    % rasters from contrast sweeps
    subplot_tight(4,4,11,tt)
    line(x{ii,7}, y{ii,7},'color','k');
    line(res.psthtimes.contpsth-res.psthtimes.contpsth(1) , p.contrastsweep.stim + p.Nrepeats+2,'color',stimcol);
    xticks(0:p.contrastsweep.duration/p.fps/4:p.contrastsweep.duration/p.fps);
    yticks(0:(p.Nrepeats+1)/2:p.Nrepeats+1);
    pbaspect([2 1 1]);           ax = gca;       ax.XColor = 'none';
    title(['Rsq contrast sweep: ',num2str(round(res.rsqall(ii,7),2))],'FontSize',9)
    
    % response per frequency for more infor check ==> Okawa, Wong, nature communication 2019, figure 7.
    subplot_tight(4,4,12,tt+0.01)
    plot(res.freqsweep.freqs, res.freqsweep.maxeachfreq(ii,:),'o','MarkerFaceColor',rgb('dodgerblue'),'Color','k')
    hold on;
    plot(res.freqsweep.fitxax,res.freqsweep.modelresp(ii,:),'color',rgb('crimson'),'LineWidth',3);
    box off;
    pbaspect([1 1 1]);     xlabel('frequencies (Hz)');     ylabel('rate');
    
    % on-off step psth
    subplot_tight(4,4,13,tt)
    line(res.psthtimes.grfirst, res.preframe.onoff(ii,:),'color',0.5* [1 1 1]);
    line(res.psthtimes.onpsth, res.onoff.onpsth(ii,:),'color',oncol);
    line(res.psthtimes.offpsth, res.onoff.offpsth(ii,:),'color',offcol);
    line(res.psthtimes.grfrq, res.preframe.freqsweep(ii,:),'color',0.5* [1 1 1]);
    xticks(round(p.changepoints(1:5)));     box off;     %  yticks(0:50:1000);
    pbaspect([2 1 1]);              xlabel('time (s)');         ylabel('rate');
    
    % frequency sweep psth
    subplot_tight(4,4,14,tt)
    line(p.freqsweep.tvec, res.freqsweep.psth(ii,:),'color',rgb('dodgerblue'));
    xticks(0:2:20);     box off;
    pbaspect([2 1 1]);              xlabel('time (s)');
    
    % contrast sweep psth
    subplot_tight(4,4,15,tt)
    line(p.contrastsweep.tvec, res.contsweep.psth(ii,:),'color',rgb('orangered'));
    xticks(0:2:20);     box off
    pbaspect([2 1 1]);              xlabel('time (s)');
    
    % contrast response curve, check ==> Okawa, Wong, nature communication 2019, figure 7
    subplot_tight(4,4,16,tt+0.01)
    plot(res.contsweep.contrasts, res.contsweep.maxeachperiod(ii,:),'o','MarkerFaceColor',rgb('orangered'),'Color','k')
    hold on;
    plot(res.contsweep.fitxax,res.contsweep.modelresp(ii,:),'color',rgb('dodgerblue'),'LineWidth',3);
    box off;
    pbaspect([1 1 1]);          xlabel('contrast');     ylabel('rate');
    
    
    if isfield(p,'sortinfo'), chinfo =  p.sortinfo(ii); else, chinfo = clus(ii,:); end
    [filename,pngfilename] = rgcname('Chirp Stimulus', chinfo, p.date, ii);
    
    suptitle(h, filename,2);
    savepngFast(h, savingpath, pngfilename);
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', ii,size(clus,1),toc);
    fprintf(msg);
    close(h);
end


end

