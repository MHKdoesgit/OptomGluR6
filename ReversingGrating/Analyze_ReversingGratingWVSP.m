

function rgwvspdata = Analyze_ReversingGratingWVSP(datapath, varargin)
%
%%% Analyze_ReversingGratingWVSP %%%
%
%
% This function analyze moving bar stimulus, the analysis includes calculation
% of the rasters and construction of DS polar plots along with DSI and OSI
% measurement. Fianlly this function calculate the p-value for significance
% of the DSI and OSI compared to a randomly shuffeled distribution of the
% spikes.
%
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   mbdata : a cell structure containing rasters values, polar plot values
%            and DSI, OSI measurement for each ganglion cell.
%   Plot : This function will also plot the rasters and polar plot for all
%          the input cells.
%
% written by Mohammad, 14.02.2018 (Valentine s day <3).


%function [res]= analyzeReversingGratingWVSPnew2(experiment, expId, varargin)
%ANALYZEREVERSINGGRATING Analyzes responses to OnOffGrating stimuli


[thisExp, savingpath] = loadRawData(datapath,'reversinggratingwvsp','ReversingGratingWVSP_Analysis');
stimPara            =       thisExp.stimPara;
stimPara.bindur     = 10/1e3;
stimPara.screen     = thisExp.info.screen.resolution;
stimPara.pixelsize  = thisExp.info.screen.pixelsize;

% if strcmpi(thisExp.lightprojection,'lightcrafter')
%     stimPara.monitorpixel = 8;
% else
%     stimPara.monitorpixel = 7.5;
% end

if ~isequal(stimPara.stripewidths(end),stimPara.screen(1)), stimPara.stripewidths(end)=stimPara.screen(1); end
stimPara.fs = double(thisExp.info.samplingrate);
%res= thisExp; 
%--------------------------------------------------------------------------
disp('Reconstructing stimulus image...');tic;
stimulus = reconstructReversingGrating(stimPara, stimPara.screen(1), stimPara.screen(2));
res.stimulus     =      stimulus;
swidths          =      stimPara.stripewidths; 
sphases          =      stimPara.Nphases;
Nwidths          =      numel(swidths); 
Nphases          =      sum(sphases);
fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
Ncells           =      size(thisExp.clusters,1);
pulsetimes       =      round(thisExp.ftimes * stimPara.fs); 
spktimes         =      double(spikeCell2Mat(thisExp.spiketimes, stimPara.fs));

phasePulses      =      (stimPara.preframes>0) + (stimPara.Nreversals)*(1+(stimPara.grayframes > 0));
Ntrials          =      floor(numel(pulsetimes)/(phasePulses*Nphases));
pulsetimes       =      pulsetimes(1:phasePulses*Nphases*Ntrials);

pulsestarts      =      reshape(pulsetimes,phasePulses,Nphases*Ntrials);
pulsestarts      =      pulsestarts((stimPara.preframes > 0)+1:end,:);
revToUse         =      floor(stimPara.Nreversals/2);
pulsestarts      =      pulsestarts(1:2*revToUse,:);
alldurs          =      diff(pulsestarts);

if stimPara.grayframes>0
    graydur     =       round(mean(mean(alldurs(1:2:end,:))));
    sdur        =       round(mean(mean(alldurs(2:2:end,:))));
    pulsestarts =       pulsestarts(2:2:end,:);
else
    graydur     =       0;
    sdur        =       round(mean(alldurs(:)));
end
trialduration   =       graydur+sdur;
%--------------------------------------------------------------------------
disp('Binning and rearranging spike times...');tic;

allRasters      =       ndSparse.build([spktimes(:,1) spktimes(:,2)],1);


%baseline calculation 
bsltime         =       0.1; 
bsldur          =       round(bsltime * stimPara.fs);
bslpulses       =       pulsestarts(1,:)-bsldur;
bslBins         =       uint32(bslpulses(:)'+(0:bsldur-1)');
bslRasters      =       allRasters(:,bslBins);
bslRates        =       reshape(bslRasters,Ncells,bsldur,Nphases*Ntrials);
bslRates        =       full(squeeze(sum(bslRates,2)))/bsltime;
baselineRates   =       mean(bslRates,2);  
baselineRatesSEM=       nansemfun(bslRates,2);
res.allBaseline =       baselineRates;

indexBins       =       uint32(pulsestarts(:)'+(0:trialduration-1)');
allRasters      =       allRasters(:,indexBins);
allRasters      =       reshape(allRasters,[Ncells,trialduration*2,revToUse,Nphases,Ntrials]);

rasterTimes     =       ((0:(2*trialduration-1))+0.5)/stimPara.fs;
res.allRasters  =       allRasters;
res.rasterTimes =       rasterTimes;
res.trialdur    =       2*trialduration;

fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
disp('Binning rasters for PSTH and calculations...');tic;

binlength       =       round(stimPara.bindur*stimPara.fs); 
psthbins        =       floor(2*trialduration/binlength);

%psths only for calculations
allRates        = reshape(allRasters(:,1:binlength*psthbins,:,:,:),...
                [Ncells, binlength, psthbins,revToUse, Nphases, Ntrials]);
allRates        = full(squeeze(sum(allRates,2)))/stimPara.bindur; %bin
allRates        = permute(allRates(:, :, 2:end, :, :), [1 3 5 2 4]); % remove first reversal
allRates        = reshape(allRates, [Ncells, Ntrials*(revToUse-1), psthbins, Nphases]);

qualityRsq      = rfroutines.imageTrialRsq(reshape(allRates,Ncells, Ntrials*(revToUse-1),[]));
res.qualityRsq  = qualityRsq;

allRatesSEM     =   squeeze(nansemfun(allRates,2)); 
res.allRatesSEM =   allRatesSEM;
allRates        =   squeeze(mean(allRates,2)); 
res.allRates    =   allRates;

allTimes        =   ((0:psthbins-1)+0.5)*stimPara.bindur; 
res.allTimes    =   allTimes; 
fprintf('Done! Took %2.2f s...\n', toc);
%--------------------------------------------------------------------------
disp('Calculating the Fourier components...');tic;

fftRates        = allRates;

L       =   size(fftRates,2);        % Length of signal
Y       =   fft(fftRates,[],2);
P2      =   abs(Y/L);
P1      =   P2(:,1:floor(L/2)+1,:);
P1(:,2:end-1,:)     =   2*P1(:,2:end-1,:);


allF1       =       NaN(Ncells,Nwidths); 
allF2       =       NaN(Ncells,Nwidths);
rateInds    =       NaN(Ncells, Nwidths,2);
allF2alt    =       NaN(Ncells, Nwidths);

for iw = 1:Nwidths
    inds    =   sum(sphases(1:iw-1))+1:sum(sphases(1:iw-1))+sphases(iw);
    
    [f1, f1ind] = max(squeeze(P1(:,2,inds)),[],2);
    allF1(:,iw) = f1;  
    rateInds(:,iw,1) = inds(f1ind);
    
    [f2, f2ind] = max(squeeze(P1(:, 3, inds)),[],2);
    allF2(:,iw) = f2; 
    rateInds(:,iw,2) = inds(f2ind);
    
    allF2alt(:,iw) = mean(P1(:, 3, inds),3);
end

res.allF1   =   allF1;
res.allF2   =   allF2; 
res.allF2mean   =   allF2alt;

allRatio        =   allF2./allF1; 
res.allRatio    =   allRatio;
res.ratioInd    =   allRatio(:,end-(swidths(end) > 100));
res.rateInds    =   rateInds;
nlinIdx         =   max(allF2,[],2)./max(allF1,[],2); 
nlinIdx(qualityRsq < 0.1)     =   NaN; 
res.nlinIdx     =   nlinIdx;
classicNlinIdx  =   max(allF2alt./allF1,[],2); 
classicNlinIdx(qualityRsq<0.1)  =   NaN;
res.classicNlinIdx      =       classicNlinIdx;

disp('Done!');toc;
%--------------------------------------------------------------------------
disp('Fitting the maximum rates...');tic;

maxRates    =   NaN(Ncells, Nwidths); 
maxRatesSEM =   NaN(Ncells,Nwidths);
% targetFrames=prestimframes+[1:stimPara.Nframes...
%     stimPara.Nframes+stimPara.grayframes+(1:stimPara.Nframes)];

for iw = 1:Nwidths
    inds    =   sum(sphases(1:iw-1))+1:sum(sphases(1:iw-1))+sphases(iw);
    trates  =   reshape(allRates(:,:,inds),Ncells,[]);
    tratessem       =   reshape(allRatesSEM(:,:,inds),Ncells,[]);
    [mrates,im]     =   max(trates,[],2);
    maxRates(:,iw)  =   mrates;
    maxRatesSEM(:,iw)   =   tratessem(sub2ind([Ncells size(tratessem,2)],(1:Ncells)',im));
end

maxRates     =  [baselineRates maxRates]; %append the baseline
res.maxRates =  maxRates; %res.maxF2=maxF2; res.meanF1=meanF1;

maxRatesSEM  =  [baselineRatesSEM maxRatesSEM]; %append the baseline
res.maxRatesSEM     =   maxRatesSEM; 

allParams    =      NaN(Ncells,4);
xx           =      [0 swidths(1:end-(swidths(end)>100))]*stimPara.pixelsize;%*1e6;

optimopts    =      optimset('Display','off','TolX', 1e3*eps,'TolFun', 1e3*eps,'MaxFunEvals',...
                              2e4,'MaxIter',2e4,'Jacobian', 'on');
for ii = 1:Ncells
    yy      =   maxRates(ii,1:end-(swidths(end)>100));
    if max(yy) < 10; continue; end
    
    %Construct the guess
    minSpikes   =   min(yy)-1;
    maxSpikes   =   max(yy)+1;
    spikesT     =   log((yy-minSpikes)./(maxSpikes-yy));
    p           =   polyfit(xx, spikesT, 1); %linear fit to get the estimates for the guess
    guess       =   [max([minSpikes 0]) maxSpikes -p(2)/p(1) max([p(1) 0])];
    try
        params  =   lsqcurvefit(@(p,x) logistic4(p,x), guess, xx,yy', [0 0 0 0], [38 1000 500 Inf], optimopts);
    catch ME
        disp(ME.message);
        params  =   lsqcurvefit(@(p,x) logistic4(p,x), guess, xx,yy');
    end
    allParams(ii,:)     =   params;
end
res.allParams   =   allParams;

%rgwvspdata.allParams = allParams;

spatialScales       =       allParams(:,3); %-log(1/options.scalePercentage-1)./allParams(:,4);
spatialScales(qualityRsq<0.1)   =    NaN;
res.spatialScales   =       spatialScales;
disp('Done!');toc;
%--------------------------------------------------------------------------
disp('Saving data...');tic;
res.para            =       stimPara;
res.savingpath      =       savingpath;
if isfield(thisExp,'sortinginfo')
    res.sortinfo    =       thisExp.sortinginfo; 
end
rgwvspdata          =       res;
filename            =       [num2str(stimPara.expnumber,'%02d'),...
    '-Reversing_grating_wvsp_analysis_for_experiment_on_',thisExp.date,'.mat'];

save([savingpath,filesep,filename],'-struct','rgwvspdata');
disp('Done!'); toc;




end