

function varargout = Analyze_LightStepsFromDarkness(datapath)
%
%%% Analyze_LightStepsFromDarkness %%%
%
%
% This function analyzes Light Steps from Darkness stimulus. This stimulus
% was specifically designed for the OptomGluR6 project. It is a copy of the
% stimulus used in the Sonja's paper. The idea is to start for black or
% zero contrast and in every step increase the brightness by 2 fold. In the
% end in 8 steps it goes from black to white. After every stimulus steps,
% the black screen is shown for same duration as the stimulus. It ain't got
% no gray screen in between.
%
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   lsfddata : a cell structure containing rasters values, psths, time-to-peak
%              decay time, etc measurements for each ganglion cell.
%   Plot : This function will also plot the rasters and psth etc for all
%          the input cells.
%
% written by Mohammad, 08.11.2018.
% updated to new version for opto project on 13.01.2021.

totaltime =tic;
if (nargin < 1),  datapath = uigetdir();       end
[lsfddata, para, clusters, savingpath] = LSFD_analysis(datapath);
% plotting
LSFD_plot(lsfddata, para, clusters, savingpath);

% Gooooooooooonnnnnggggggggggggggg!!!!
sound(struct2array(load('gong.mat','y')))
disp(seconds2human (toc(totaltime)));
varargout{1} = lsfddata;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [lsfddata, para, clus, savingpath] = LSFD_analysis(datapath,varargin)

[thisExp, savingpath] = loadRawData(datapath,'lightstepsfromdarkness','LightStepsFromDarkness_Analysis');


clus = thisExp.clusters;

para = rmfield(thisExp.stimPara,{'lmargin','rmargin','tmargin','bmargin'});
para.binlength = 10/1e3;
conts = 1;
for ii = 1:para.Nintensities-1, conts = [conts, conts(end)/2]; end  %#ok there should be a better way, fix later!
para.contraststeps = conts;
para.fps = thisExp.info.screen.refreshrate;
stimdur = para.Nframes / para.fps;
trialdur = stimdur*2;
pfrdur = para.preframes / para.fps;
startstim = pfrdur/2:trialdur:pfrdur+trialdur+para.Nsteps;
if numel(startstim) > para.Nintensities,  startstim = startstim(1:para.Nintensities);   end
para.stimonset = startstim;

ft = reshape(thisExp.ftimes(1:(1+para.Nsteps*2)*para.Nintensities*para.Nrepeats ), (1+para.Nsteps*2),[]);
% pulsetimes = reshape(thisExp.ftimes(1:(1+para.nreversals)*para.numphases*para.numtrials ), (1+para.nreversals),[]);

bkgtimes=ft(2,:); % background comes only between the intensities on the repeats or Nsteps! (shitty stimulus design)
ft=ft(2:2:end,:);
ft=ft(1:end,:);
trialdur = mean(mean(diff(ft)));

onsettimes = reshape(ft,size(ft,1),para.Nintensities,[]);
onsettimes = permute(onsettimes,[1 3 2]);
onsettimes = reshape(onsettimes,size(ft,1)*para.Nrepeats,[])';
%--------------------------------------------------------------------------
disp('Organizing the rasters...');
tic;
allRasters = createTrialRasters(thisExp.spiketimes,onsettimes,onsettimes+trialdur);

rasconts = cell(size(clus,1),para.Nintensities);
rasters = cell(size(clus,1),1);
for ii = 1:size(clus,1)
    thisrasters = [];
    tt = para.stimonset(1);
    for jj = 1: para.Nintensities
        rasconts{ii,jj} = CelltoMatUE(squeeze(allRasters(ii,jj,:)));
        thisrasters = [thisrasters, tt+ rasconts{ii,jj}];   %#ok
        tt = tt+trialdur;
    end
    rasters{ii} = thisrasters;
end
lsfddata.rastersall = rasters;
lsfddata.rasters = rasconts;
toc;
%--------------------------------------------------------------------------
disp('Estimating quality Rsq...');
tic;
ratesOdd = calculatePSTH (allRasters(:,:,1:2:end),[0 trialdur], 2*para.binlength);
ratesEven = calculatePSTH (allRasters(:,:,2:2:end),[0 trialdur], 2*para.binlength);
ratesOdd = reshape(ratesOdd, size(clus,1),[]);
ratesEven = reshape(ratesEven, size(clus,1),[]);
qualityRsq = 1-sum((ratesOdd-ratesEven).^2,2)./sum((ratesOdd-repmat(mean(ratesOdd,2),[1 size(ratesOdd,2)])).^2,2);
lsfddata.qualityRsq = qualityRsq;
toc;
%--------------------------------------------------------------------------
%because of the raster nature, easy to apply PSTH to all
[allRates,allTimes] = calculatePSTH(reshape(allRasters, size(clus,1)* para.Nintensities,[]), [0 trialdur], para.binlength, 0);
allRates = reshape(allRates,size(clus,1), para.Nintensities,[]);
lsfddata.psths = allRates(:,:,1:floor(size(allRates,3)/2)*2);
lsfddata.psthtime = allTimes(1:floor(size(allRates,3)/2)*2); %res.allSigmas=allSigmas;
%--------------------------------------------------------------------------
disp('Estimating baseline activity...');
tic;
pfrdur = (para.preframes/para.fps);
bkgRasters = createTrialRasters(thisExp.spiketimes,bkgtimes-(trialdur*2),bkgtimes);
bkgrasall = cell(size(clus,1),1);
for ii = 1:size(clus,1)
    bkgrasall{ii,1} = CelltoMatUE(bkgRasters(ii,:));
end
[bkgRates, bkgtimes] = calculatePSTH (bkgRasters,[0 pfrdur], para.binlength);
lsfddata.bkgrasters = bkgrasall;
lsfddata.bkgpsth = bkgRates(:,1:end);
lsfddata.bkgpsthtime = bkgtimes;
bsl = cellfun(@numel,bkgRasters)/(pfrdur);
bsl = mean(bsl,2);
lsfddata.baselinerate = bsl;
toc;
%---------------------------------------------------------------------------
disp('Calcualting peak latency, peak value and decay time...');
tic;
[onlat, offlat, onpeaks, offpeaks, ondecaytime, offdecaytime] = deal(nan(size(clus,1),para.Nintensities));
[ondecayprop, offdecayprop] = deal(cell(size(clus,1),2));

warning('off');
for ii = 1:size(clus,1)
    thispsth = squeeze(lsfddata.psths(ii,:,:))';
    onpsth = thispsth(1:size(thispsth,1)/2,:);
    offpsth = thispsth(1+size(thispsth,1)/2:end,:);
    
    [onlat(ii,:), onpeaks(ii,:)] = psthPeakLatency(onpsth,[25 500],10, stimdur*1e3);
    [offlat(ii,:), offpeaks(ii,:)] = psthPeakLatency(offpsth,[25 500],10, stimdur*1e3);
    
    [ondecaytime(ii,:), ondcval,ondctvec] = psthDecayTime(onpsth, onlat(ii,:), stimdur*1e3);
    
    [offdecaytime(ii,:), offdcval,offdctvec] = psthDecayTime(offpsth, offlat(ii,:), stimdur*1e3);
    
    ondecayprop{ii,1} = CelltoMatUE(ondctvec);
    ondecayprop{ii,2} = CelltoMatUE(ondcval);
    
    offdecayprop{ii,1} = CelltoMatUE(offdctvec);
    offdecayprop{ii,2} = CelltoMatUE(offdcval);
    
end
warning('on');
para.date = thisExp.date;
lsfddata.onpeaklatency = onlat;
lsfddata.onpeak = onpeaks;
lsfddata.offpeaklatency = offlat;
lsfddata.offpeak = offpeaks;
lsfddata.ondecaytime = ondecaytime;
lsfddata.ondecayprop = ondecayprop;
lsfddata.offdecaytime = offdecaytime;
lsfddata.offdecayprop = offdecayprop;
lsfddata.clusters = clus;
lsfddata.para = para;

% saving data
filename = [num2str(para.expnumber,'%02d'),'-LightStepsFromDarkness_analysis_for_experiment_on_',para.date,'.mat'];
save([savingpath,filename],'-struct','lsfddata');

toc;

end

%--------------------------------------------------------------------------------------------------%

function LSFD_plot(lsfd, para, clus, savingpath)

stimdur = para.Nframes / para.fps;
trialdur = stimdur*2;
pfrdur = para.preframes / para.fps;
startstim = pfrdur/2:trialdur:pfrdur+trialdur+para.Nsteps;
if numel(startstim) > para.Nintensities,  startstim = startstim(1:para.Nintensities);   end
% plotting modules
patfun = @(n1,n2,n3,n4,c)(patch([repmat( n1,1,2) repmat(n2,1,2)],[n3, n4, n4, n3],0,'facecolor',c,...
    'edgecolor','k','linewidth',0.25,'linestyle','-'));
pltfun = @(x,y,c1,c2,varargin)(plot(x,y,'-o','color',rgb(c1),'markerfacecolor',rgb(c2),...
    'linewidth',2,'markersize',7,varargin{:}));
xAx = trialdur*para.Nintensities+stimdur*2;
cols = lines(para.Nintensities);
msg=[];

for ii = 1:size(clus,1)
    % start plotting
    h = figure('pos',[50 50 1500 950],'color','w','vis','off');
    
    % rasters
    subplot_tight(3,6,[1 5],0.03)
    rasterPlotter(lsfd.rastersall{ii}',[],'k');
    ntrials = size(lsfd.rastersall{ii},1);
    recheight = ceil(ntrials*0.15); % para.Nsteps*para.Nrepeats+1
    
    axis([0 xAx 0 ntrials]);
    %rectangle('pos', [0 ntrials+2 xAx recheight-2],'facecolor','k','edgecolor','none');
    for kk =1:para.Nintensities
        patfun(startstim(kk),startstim(kk)+stimdur,ntrials+2,ntrials+recheight,para.contraststeps(1+para.Nintensities-kk)*[1,1,1]);
        line([startstim(kk) startstim(kk)],[0.4 ntrials+0.4],'color','r','linestyle','--');     
        %xline(startstim(kk)+stimdur,'--r');
    end
    axis([startstim(1) xAx 0 ntrials+recheight+2]);
    ax = gca;                           ax.XColor = 'none';
    ax.XTick = [];   ylabel('Trials');      ax.TickLength = [0.0025 0.0025];
    title(['Rasters ( ',num2str(para.Nsteps),' by ',num2str(para.Nrepeats),' trials for ',num2str(para.Nintensities),...
        ' intensities, ( Rsq: ',num2str(round(lsfd.qualityRsq(ii),2)),' ))']);
    ax.Position = [ax.Position(1) ax.Position(2)+0.02 ax.Position(3:4)];
    set(gca,'fontsize',8);      ax.YTick = 0:ntrials/2:ntrials;
    
    % psth
    subplot_tight(3,6,[7 11],0.03)
    p = squeeze(lsfd.psths(ii,:,:))';
    yAx = ceil(max(p(:))/10)*10;    if yAx==0, yAx=10; end
    recheight = ceil(yAx*0.1);
    for kk =1:para.Nintensities
        plot(startstim(kk)+lsfd.psthtime,p(:,kk),'color','k');
        hold on;
        patfun(startstim(kk),startstim(kk)+stimdur,yAx+2,yAx+recheight+2,para.contraststeps(1+para.Nintensities-kk)*[1,1,1]);
        line([startstim(kk) startstim(kk)],[0 yAx],'color','r','linestyle','--');
    end
    box off;
    plot(lsfd.onpeaklatency(ii,:)/1e3+startstim,lsfd.onpeak(ii,:),'mx','linewidth',2,'MarkerSize',10);
    plot(lsfd.offpeaklatency(ii,:)/1e3+startstim+stimdur,lsfd.offpeak(ii,:),'bo','linewidth',2,'MarkerSize',7);
    axis([startstim(1) xAx 0 yAx+recheight+2]);
    ax = gca;        ax.TickLength = [0.0025 0.0025];
    %rectangle('pos',[0 yAx+2 xAx recheight],'facecolor','k','edgecolor','none');
    xlabel('Time (sec)');       ylabel('Firing rate (Hz)');
    ax.Position = [ax.Position(1) ax.Position(2)+0.04 ax.Position(3:4)];
    set(gca,'fontsize',8);      ax.YTick = 0:yAx/2:yAx;
    
    % preframe rasters
    subplot_tight(3,6,6,0.03)
    rasterPlotter(lsfd.bkgrasters{ii},[],cols(1,:));
    tt = 0;
    m = mod(1:para.Nintensities * para.Nrepeats,para.Nintensities); m(m==0)= para.Nintensities;
    for qq = 1: para.Nintensities * para.Nrepeats
        rectangle('pos',[pfrdur+0.2 0.6+tt 0.5 0.8],'facecolor',para.contraststeps(1+para.Nintensities-m(qq))*[1,1,1]);
        tt = tt+1;
    end
    axis([0 pfrdur+0.71 0 para.Nintensities * para.Nrepeats+1]);
    title('Preframe rasters');  box off;    axis square;
    yticks(0:para.Nintensities:para.Nintensities*para.Nrepeats);    xlabel('Time (sec)');
    %ax = gca;   ax.Position = [ax.Position(1) ax.Position(2)-0.04 ax.Position(3:4)];    
    set(gca,'fontsize',8);
    
    % preframe PSTH
    subplot_tight(3,6,12,0.03)
    area(lsfd.bkgpsthtime,lsfd.bkgpsth(ii,:),'facecolor',rgb('silver'),'EdgeColor',rgb('dimgray'),'LineWidth',0.01);
    yline(lsfd.baselinerate(ii,:),'--r','baseline');
    xlabel('Time (sec)');       box off;        axis square;        title('Preframe PSTH');
    yAxpfr = ceil(max(lsfd.bkgpsth(ii,:))/2)*2;    if yAxpfr==0, yAxpfr=10; end
    axis([0 pfrdur 0 yAxpfr]);     yticks(0:yAxpfr/2:yAxpfr);        xticks(0:1:20);
    ax = gca;   ax.Position = [ax.Position(1) ax.Position(2)+0.01 ax.Position(3:4)];    set(gca,'fontsize',8);
    
    % peak values
    subplot_tight(3,3,7,0.06)
    pltfun(1:para.Nintensities, lsfd.offpeak(ii,:),'royalblue','deepskyblue');
    hold on;
    pltfun(1:para.Nintensities, lsfd.onpeak(ii,:),'crimson','red');
    yline(lsfd.baselinerate(ii,:),'--k','baseline');
    axis square;    axis([0.5 para.Nintensities+0.5 0 yAx]);        box off;
    xticklabels(round(fliplr(para.contraststeps*100),2));   xticks(1:1:para.Nintensities);
    yticks(0:yAx/2:yAx);    xlabel('Contrast values');  ylabel('Firing rate (Hz)');
    title('Peak values');   legend('OFF steps','ON steps','location','northeastoutside');  legend boxoff;
    set(gca,'fontsize',8); pbaspect([1.5 1 1]);
    
    % time-to-peaks
    subplot_tight(3,3,8,0.06)
    pltfun(1:para.Nintensities, lsfd.offpeaklatency(ii,:),'cyan 3','cyan')
    hold on
    pltfun(1:para.Nintensities, lsfd.onpeaklatency(ii,:),'salmon','pink');
    axis([0.5 para.Nintensities+0.5 0 250]);    axis square;        box off;
    xticklabels(round(fliplr(para.contraststeps*100),2));   xticks(1:1:para.Nintensities);
    yticks(0:50:1000);    xlabel('Contrast values');  ylabel('Time-to-peak (ms)');
    title('Time-to-peak');   legend('OFF steps','ON steps','location','northeastoutside'); legend boxoff;
    set(gca,'fontsize',8);  pbaspect([1.5 1 1]);
    
    % decay time
    subplot_tight(3,3,9,0.06)
    pltfun(1:para.Nintensities, lsfd.offdecaytime(ii,:),'teal','lightseagreen');
    hold on;
    pltfun(1:para.Nintensities, lsfd.ondecaytime(ii,:),'peru','melon');
    axis square;    axis([0.5 para.Nintensities+0.5 0 250]);
    xticklabels(round(fliplr(para.contraststeps*100),2));   xticks(1:1:para.Nintensities);
    yticks(0:50:1000);    xlabel('Contrast values');  ylabel('Decay time (ms)');
    title('Decay time');   legend('OFF steps','ON steps','location','northeastoutside'); legend boxoff;
    set(gca,'fontsize',8); box off; pbaspect([1.5 1 1]);
    
    % saving plot
    if isfield(p,'sortinfo')
        [filename,filenamesavepng] = rgcname('Light Step from Darkness', para.sortinfo(ii),savingpath);
        pngfilename = [num2str(ii,'%02d-'),extractAfter(filenamesavepng,[savingpath,'\'])];
    else
        filename = rgcname('Light Step from Darkness', clus(ii,:),savingpath);
        pngfilename = [num2str(ii,'%02d-'),filename];
    end
    suptitle(h,filename,3);
 
    savepngFast(h,savingpath,pngfilename);
    close(h);
    
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', clus(ii,1),clus(ii,2));
    fprintf(msg);
    
end

end
