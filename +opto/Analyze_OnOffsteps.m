

function varargout = Analyze_OnOffsteps(datapath, varargin)
%
%%% Analyze_OnOffsteps %%%
%
%
% This function measure the spikes during each color stimulation and plot
% the rasters and psth for each cell. Its an update to the previos RGBsteps
% and RGBPlotter.
%
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder (include \Data_Analysis\).
%
%================================Output====================================
%
%   RGBData : a cell structure containing stimulus parameters, for
%              four color stimuluation conditions and thier coresponding
%              PSTH also it measure the reponse type and cell type of
%              each recoded ganglion cell.
%   Plot : This function will also plot raters and psth plot for all 4
%          chromatic stimulus conditions.
%
%
% written by Mohammad, 15.03.2015
% fixed the bug for histc and change scale to sec on 12.01.2017.
% updated to completely new version with better plot design on 08.02.2018.
% updated to new version to handel iPRGC stim and have better folder naming
% on 09.05.2019.
% updated to new version with single file saved and better plotting.

totaltime =tic;
if (nargin < 1),    datapath = uigetdir();  end

thisExp = loadRawData(datapath,{'rgbsteps','onoffsteps'});

%%%folder making
stimname = lower(thisExp.stimPara.originalname);
if contains(stimname,'rgbsteps')
    savefolder = 'RGBsteps_Analysis';
    %dsfold = 'rgb_data';
elseif contains(stimname,'onoffsteps') && ~contains(stimname,'iprgc')
    savefolder = 'OnOffsteps_Analysis';
    %dsfold = 'onoff_data';
elseif contains(stimname,'iprgc')
    savefolder = 'iPRGC_Analysis';
    %dsfold = 'iprgc_data';
end

savingPath = [datapath,'/Data Analysis/', sprintf('%02G', thisExp.stimPara.expnumber),'-',savefolder,'/'];
if ~exist(savingPath,'dir'), mkdir(savingPath); end
%if ~exist([savingPath,dsfold],'dir'), mkdir ([savingPath,dsfold]);  end
clearvars -except thisExp totaltime savingPath datapath dsfold savefolder;

para = thisExp.stimPara;
expinfo = thisExp.info;
if isfield(para,'usered')
    para = rmfield(para,{'lmargin','rmargin','tmargin','bmargin','usered','usegreen',...
        'useblue','redmeanintensity','greenmeanintensity','bluemeanintensity'});
else
    para = rmfield(para,{'lmargin','rmargin','tmargin','bmargin'});
end
% analysis variables
if strcmpi( savefolder, 'iPRGC_Analysis')
    para.binlength = 200/1e3;
else
    para.binlength = 10/1e3;
end
if (thisExp.stimPara.preframes > 0), ftgap = 4; else, ftgap = 2; end
para.fps = expinfo.screen.refreshrate;
para.nbins = (thisExp.stimPara.Nframes+thisExp.stimPara.preframes)*2/para.fps*(1/para.binlength);
if isfield(thisExp,'sortinginfo')
    para.sortinfo = thisExp.sortinginfo;
end
para.date = thisExp.date;


rasall = cell(size(thisExp.clusters,1),1);
psthall = zeros(size(thisExp.clusters,1),para.nbins);
msg = [];

for ii = 1: size(thisExp.clusters,1)
    [rasters, psth] = get_rasters_psth(thisExp.ftimes,thisExp.spiketimes{ii}, ftgap, para.binlength);
    if isfield(thisExp,'sortinginfo'), si = para.sortinfo(ii); else, si = []; end
    tvec = plot_onoff_steps(rasters,psth, para, thisExp.clusters(ii,:), savingPath, si, ii); %para.sortinfo(ii)
    rasall{ii} = rasters;
    psthall(ii,:) = psth;
    %paraall{ii} = rgbdat.para;
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', ...
        thisExp.clusters(ii,1),thisExp.clusters(ii,2));
    fprintf(msg);
end
para.tvec = tvec;

onoff.rasters = rasall;
onoff.psth = psthall;
if isfield(thisExp,'sortinginfo')
    onoff.sortinfo = para.sortinfo;
    onoff.para = rmfield(para,'sortinfo');
else
    onoff.para = para;
end
filename = [upper(savefolder(1)),lower(savefolder(2:end)),' for experiment on ', para.date];
filename = [num2str(para.expnumber,'%02d'),strrep(filename,'_',' ')];
save([savingPath,'\',filename,'.mat'],'-v7.3','-struct','onoff');
sound(struct2array(load('gong.mat','y')));
disp(seconds2human (toc(totaltime)));
varargout{1} = onoff;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [rasters, psth] = get_rasters_psth(ft,spk, steps, binLength, varargin)

stim = ft(1:steps:numel(ft));
stimround = round(mean(diff(stim))); % one --|'|--|_| is 5000 ms, 4 frame change
binVec = linspace(0,stimround,1+stimround/binLength);
rasTrials = cell((length(stim)-1),1);

for i = 1:numel(stim)-1 %first trial is background intensity
    rasTrials{i} = spk(and(spk > stim(i),spk <= stim(i+1))) - stim(i);
    if isempty(rasTrials{i}); rasTrials{i} = NaN; end
end
rasters = CelltoMatUE(rasTrials);
if size(rasters,2) < 2
    psth = (histc(rasters',binVec)/size(rasters,1))*(1/binLength);  % to fix situation with 1 spike
else
    psth = mean(histc(rasters',binVec),2)*(1/binLength);   % must transpose contspk for histc to work properly.
end
psth = psth(1:end-1);

end

%--------------------------------------------------------------------------------------------------%

function varargout = plot_onoff_steps(rasters, psth, para, clus, savingPath, sortinfo, iter, varargin)

pspltt0 = @(t1,t2,p,c)(plot(linspace(t1,t2,diff([t1,t2])/para.binlength)-p, ...
    psth(int16(1+t1/para.binlength: t2/para.binlength)),'color',c,'linewidth',2,varargin{:}));
psplt = @(t1,t2,p,c)(plot(linspace(t1,t2,diff([t1,t2])/para.binlength+1)-p, ...
    psth(int16(t1/para.binlength: t2/para.binlength)),'color',c,'linewidth',2,varargin{:}));


stimDur = (para.Nframes/para.fps);
prfDur = (para.preframes/para.fps);
%[colorSet,stimname,titr] = colorSelector(para.originalname);
tvec = [0,prfDur,stimDur+prfDur,2*prfDur+stimDur,2*prfDur+2*stimDur];

[pfron,onspk] = early_late_spks(rasters,[tvec(1),tvec(2)],[tvec(2),tvec(3)]);
[pfroff,offspk] = early_late_spks(rasters,[tvec(3),tvec(4)],[tvec(4),tvec(5)]);
pfronneg250 = early_late_spks(rasters,[tvec(2)-0.25,tvec(2)],[tvec(2),tvec(3)]);
pfroffneg250 = early_late_spks(rasters,[tvec(4)-0.25,tvec(4)],[tvec(4),tvec(5)]);
[~,pfronpos250] = early_late_spks(rasters,[tvec(2),tvec(3)],[tvec(3),tvec(3)+0.25]);
[~,pfroffpos250] = early_late_spks(rasters,[tvec(4),tvec(5)],[tvec(1),tvec(1)+0.25]);
% plotting
oncol = rgb('crimson');
offcol = rgb('dodgerblue');
pfrcol = rgb('silver');

h = figure('pos',[20 0 1800 1050],'color','w','vis','off');
ntr = ceil(size(rasters,1)/2)*2;
yAx = max([1,ceil(max(psth)/2)*2]);
xgap = 0.05;
subplot_tight(4,4,[1,5],[0.05,0.025])
rasterPlotter(onspk',[],oncol);
if para.preframes>0
    rasterPlotter(pfronneg250',[],pfrcol);
    rasterPlotter(pfronpos250',[],pfrcol);
end
axis([tvec(2)-0.25-xgap,tvec(3)+0.25+xgap 0 ntr]);
set(gca,'xtick',tvec(2):stimDur/2:tvec(3),'ytick',0:ntr/2:ntr,'fontsize',9);
xlabel('Time (sec)'); ylabel('Trials'); title('Rasters On stimulus'); axis square;

subplot_tight(4,4,[2,6],[0.05,0.025])
rasterPlotter(offspk',[],offcol);
if para.preframes>0
    rasterPlotter(pfroffneg250',[],pfrcol);
    rasterPlotter(tvec(end)+pfroffpos250',[],pfrcol);
end
axis([tvec(4)-0.25-xgap,tvec(5)+0.25+xgap 0 ntr]);
set(gca,'xtick',tvec(4):stimDur/2:tvec(5),'ytick',0:ntr/2:ntr,'ycolor','w','fontsize',9);
xlabel('Time (sec)'); title('Rasters Off stimulus'); axis square;

subplot_tight(4,4,[9,12],0.035)
rasterPlotter(onspk',[],oncol);
rasterPlotter(offspk',[],offcol);
if para.preframes>0
    rasterPlotter(pfron',[],pfrcol);
    rasterPlotter(pfroff',[],pfrcol);
end
rectangle('pos',[tvec(2) -5 tvec(3)-tvec(2) 2],'facecolor',rgb('gold'),'edgecolor','none');
rectangle('pos',[tvec(4) -5 tvec(5)-tvec(4) 2],'facecolor',rgb('black'),'edgecolor','none');
if para.preframes>0
    rectangle('pos',[tvec(1) -5 tvec(2) 2],'facecolor',pfrcol,'edgecolor','none');
    rectangle('pos',[tvec(3) -5 tvec(4)-tvec(3) 2],'facecolor',pfrcol,'edgecolor','none');
end
axis([tvec(1) tvec(end) -5 ntr]);        axis off;

subplot_tight(4,4,[3,8],0.07)
hold on; box off;
psplt(tvec(2),tvec(3),0,oncol);
psplt(tvec(4),tvec(5),tvec(3),offcol);
if para.preframes>0
    psplt(tvec(2)-0.25,tvec(2),0,rgb('cadmiumorange'));
    psplt(tvec(3),tvec(3)+0.25,0,rgb('cadmiumorange'));
    psplt(tvec(4)-0.25,tvec(4),tvec(3),rgb('deepskyblue'));
    pspltt0(tvec(1),tvec(1)+0.25,-tvec(3),rgb('deepskyblue'));
    ShadePlotForEmpahsis({[tvec(2)-0.25-xgap,tvec(2)],[tvec(3),tvec(3)+0.25+xgap]},{pfrcol,pfrcol},{0.3,0.3});
end
axis([tvec(2)-0.25-xgap,tvec(3)+0.25+xgap 0 yAx]);
set(gca,'xtick',tvec(2):stimDur/2:tvec(3),'ytick',0:yAx/2:yAx,'xticklabel',0:stimDur/2:stimDur,'fontsize',9);
xlabel('Time (sec)'); ylabel('Spike rate (Hz)'); title('PSTH On & Off stimulus');

subplot_tight(4,4,[13,16],[0.045,0.035])
hold on;
p1 = psplt(tvec(2),tvec(3),0,oncol);
p2 = psplt(tvec(4),tvec(5),0,offcol);
if para.preframes>0
    pspltt0(tvec(1),tvec(2),0,pfrcol);
    psplt(tvec(3),tvec(4),0,pfrcol);
end
box off;
axis([tvec(1),tvec(end) 0 yAx]);
set(gca,'xtick',unique(tvec),'ytick',0:yAx/2:yAx,'ticklength',[0.0025 0.0025],'fontsize',9);
xlabel('Time (sec)'); ylabel('Spike rate (Hz)');
legend([p1,p2],'On','Off'); legend boxoff;

if ~isempty(sortinfo), chinfo = sortinfo; else, chinfo = clus; end
[filename,pngfilename] = rgcname([upper(para.stimulus(1)),para.stimulus(2:end)], chinfo, para.date, iter);

suptitle(h,filename,2);
savepngFast(h,savingPath,pngfilename);
close(h);

varargout{1} = tvec;
varargout{2} = filename;

end

%--------------------------------------------------------------------------------------------------%

function varargout = early_late_spks(spk,early,late,varargin)

earlyspk = spk;
latespk = spk;
earlyspk(not(spk > early(1) & spk <= early(2))) = NaN;
latespk(not(spk > late(1) & spk <= late(2))) = NaN;
varargout{1} = earlyspk;
varargout{2} = latespk;

end
