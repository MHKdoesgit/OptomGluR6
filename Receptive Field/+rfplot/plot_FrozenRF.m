

function rfdata = plot_FrozenRF(datapath, varargin)
%
%%% plot_FrozenRF %%%
%
%
% This function plot the analyzed data for the checkerflicker frozen noise
% stimulus.
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder. It can also be the analyzed data
%              struture from the Analyze_FrozenCheckerFlicker function.
%
%================================Output====================================
%
%   rfdata : output is the structure of the analyzed data.
%   Plot : Receptive field, temporal components, spatial components, Auto-
%          correlogram, surround index, response perdication and more plots.
%
% writetn by Mohammad for the opto project 11.01.2021.

totaltime = tic;

if isstruct(datapath)
    rfdata = datapath;
else
    rfdp = dir([datapath,filesep,'Data Analysis',filesep,'*FrozenNoise*']);
    
    if size(rfdp,1) > 1
        selectmsg = 'Too many frozen noise stimuli on the dance floor!! select what you want to analyze!! ==> ';
        for ii=1:   size(rfdp,1), disp(rfdp(ii).name);        end
        expnum = input(selectmsg);
        rfdp = rfdp(expnum);
    end
    
    rfsavingpath = [rfdp.folder,filesep,rfdp.name,filesep];
    rfdatpath = dir([rfsavingpath,'*checkerflicker_analysis*']);
    
    if size(rfdatpath,1) > 1
        selectmsg = 'Too many frozen noise data on the dance floor!! select what you want to analyze!! ==> ';
        for ii=1:   size(rfdatpath,1), disp(rfdatpath(ii).name);        end
        expnum = input(selectmsg);
        rfdatpath = rfdatpath(expnum);
    end
    rfdata = load([rfdatpath.folder,filesep,rfdatpath.name]);
end

rd      = rfdata;
nrows   = 4;
ncols   = 6;
tt      = 0.035;
col     = rgb('crimson');
lcol    = lines(10);
msg     = [];

for ii = 1: size(rd.clusters,1)
    
    h = figure('pos',[100 0 1750 1000],'color','w', 'vis','off');
    
    % sta best frame for green----------------------------------------------------------------------
    subplot_tight(nrows, ncols, [1 8], tt+0.02)
    sta = squeeze(rd.sta(ii,:,:,:));
    [~,mx] = max(squeeze(rd.modeltcomps(ii,:)));
    [~,mn] = min(squeeze(rd.modeltcomps(ii,:)));
    if mx > mn, peakpos = mx; else, peakpos = mn; end
    rfplot.bestframe(sta,[0 0 peakpos], squeeze(rd.contourpoints(ii,:,:)), squeeze(rd.ellipsepoints(ii,:,:)),...
        rd.moransI(ii,1),rd.para,{'color','k','linecolor',lcol(3,:),'colorbarpos',[0.04 0.012]});
    axis tight off
    t = title(['rf diameter: ',num2str(round(rd.rfdiameters(ii,1)*1e6,2)),' (µm), contour area: ',...
        num2str(round(rd.contourareas(ii,1)*1e4,2)),' (mm^2)'],'fontsize',9);
    t.Position(2) = t.Position(2)-40;
    
    % temporal components------------------------------------------------------------------
    subplot_tight(nrows,ncols,[3 4],tt)
    tc = rd.temporalComps(ii,:)./ sqrt(sum(rd.temporalComps(ii,:).^2,2)); % eucleadian normalization
    tcm = rd.modeltcomps(ii,:)./ sqrt(sum(rd.modeltcomps(ii,:).^2,2)); % eucleadian normalization
    
    plot(rd.tcTimeVec,tcm,'color','k','linewidth',0.5); hold on;
    plot(rd.tcTimeVec,tc,'color',col,'linewidth',1.5);
    xline(rd.tcTimeVec(peakpos),'--',[num2str(-rd.tcTimeVec(peakpos)*1e3),' ms'],'color',lcol(3,:));
    xlim([-0.5 0]);            box off;          pbaspect([2 1 1]);
    xticks(-1:0.2:0);     yticks(-2:0.2:2);     ylabel('filter');  %title('temporal components');
    legend('tempcomp-fit','sta (best pixel)','Location','northwest'); legend boxoff;
    
    % spatial components ---------------------------------------------------------------------------
    subplot_tight(nrows,ncols,10,tt)
    sp = squeeze(rd.spatialComps(ii,rd.rangey{ii},rd.rangex{ii}));
    mx = max(abs(sp),[],'all');  if isempty(mx), mx = 1; end;  if isnan(mx), mx = 1; end
    xsx = rd.spaceVecX(rd.rangex{ii});      xsx = xsx(~isnan(xsx));
    ysy = rd.spaceVecY(rd.rangey{ii});      ysy = ysy(~isnan(ysy));
    if length(xsx) > size(sp,2), xsx = xsx(1:size(sp,2));end
    if length(ysy) > size(sp,1), ysy = ysy(1:size(sp,1));end
    imagesc(xsx,ysy,sp);                hold on;
    rfplot.coordinates(squeeze(rd.contourpoints(ii,:,:)),'allcolors','k','lineslinewidth',1,'circlelinewidth',2)
    axis on;        xticks([]); yticks([]);
    if tcm(peakpos) < 0
        colormap(gca,cbrewer('div','RdBu',255));
    else
        colormap(gca,flipud(cbrewer('div','RdBu',255)));
    end
    caxis([-mx mx]);
    axis tight equal;   % title('spatial component');
    
    % nonlinearities--------------------------------------------------------------------------------
    subplot_tight(nrows,ncols,9,tt)
    plot(rd.nlxmodel(ii,:),rd.nlymodel(ii,:),'.-','color','k','LineWidth',1); hold on;
    plot(rd.nlx(ii,:),rd.nly(ii,:),'color',col,'LineWidth',2);
    xlim([-3 3]);       xticks(-5:1:5);       box off;    axis square;      ylabel('rate (Hz)');
    %title('nonlineaities');
    
    % all receptive fields zoom---------------------------------------------------------------------
    subplot_tight(nrows,ncols,5,tt)
    cp = permute(squeeze(rd.contourpoints(:,:,:)),[2 3 1]);
    rfplot.allRFs(cp,ii,rd.para.screen);     axis tight equal;
    xticks(0:100:rd.para.screen(1));     yticks(0:100:rd.para.screen(2));
    
    % all receptive fields -------------------------------------------------------------------------
    subplot_tight(nrows,ncols,6,tt)
    rfplot.allRFs(cp,ii,rd.para.screen,'outline',true);     axis equal;
    axis([0 rd.para.screen(1) 0 rd.para.screen(2)]);
    
    % auto-correlogram------------------------------------------------------------------------------
    subplot_tight(nrows,ncols,11,tt)
    area(rd.autoCorrLag,rd.autoCorrelations(ii,:),'facecolor',rgb('deepskyblue'),'edgecolor','none')
    xline(2,'k--');         % pbaspect([2 1 1]);
    xlim([0 50]);      xticks(0:25:100);         yticks(0:0.01:1);      box off;
    title('Auto-correlogram');
    % xlabel('time (ms)');
    
    % surround index--------------------------------------------------------------------------------
    subplot_tight(nrows,ncols,12,tt)
    plot(rd.surrSigmaVals, rd.surrActivation(ii,:),'LineWidth',1);
    axis([0 8 0 1]);    xticks(0:2:10);     yticks(0:0.5:1);    box off;     %axis square;
    %xlabel('rf sigma');
    title(['Surround index: ' ,num2str(round(rd.surroundIndex(ii,1),2))]);
    
    % sta frames -----------------------------------------------------------------------------------
    subplot_tight(nrows,ncols,[13 17],[0.01 0.01])
    rfplot.staframes(sta,20,2,rd.para,'stapeak',[0 0 peakpos],'gap',2, 'peakframeoutlinecolor',...
        col,'fontcolor', 'k','outlinecolor',rgb('dimgray'),'peakframefontcolor',col);
    ylabel('Sta frames');
    
    % trial rasters---------------------------------------------------------------------------------
    subplot_tight(nrows,ncols,18,tt)
    imagesc(squeeze(rd.frozenTrialRates(ii,:,:))');
    colormap(gca,flipud(gray));
    axis square off;        title('Rasters');
    caxis([0 max(squeeze(rd.frozenTrialRates(ii,:,:)),[],'all')/2+1]);
    
    % prediction for 2-color combined nonlinearities------------------------------------------------
    subplot_tight(nrows,ncols,[19 24],tt)
    tvec = rd.frozenTimeVec(rd.para.nonlinBinN:end);
    tf = @(r1,r2,c1,c2)(title(['Rsq: ',num2str(round(r1,2)),' ,model Rsq: ',num2str(round(r2,2)),...
        ', & cc norm: ',num2str(round(c1,2)),' , model cc norm: ', num2str(round(c2,2))]));
    
    plot(tvec,rd.frozenRates(ii,:),'k','LineWidth',0.5);          hold on;
    plot(tvec,rd.predictions(ii,:),'color',lcol(1,:),'LineWidth',1);
    plot(tvec,rd.modelpredictions(ii,:),'-','color',lcol(3,:),'LineWidth',1)
    xlim([tvec(1) tvec(end)]);      box off;    xticks(0:5:50);
    set(gca,'TickLength',[0.0025 0.0025]);
    legend('cell response','LN-model (tempcomp)','LN-model (fit tempcomp)','NumColumns',3); legend boxoff;
    tf(rd.predRsq(ii),rd.modelpredRsq(ii),rd.predCCnorm(ii),rd.modelpredCCnorm(ii));
    ylabel('Rate (Hz)');
    
    % title
    if isfield(rd.para,'sortinfo'), chinfo =  rd.para.sortinfo(ii); else, chinfo = rd.clusters(ii,:); end
    [filename,pngfilename] = rgcname('Receptive Field', chinfo, rd.para.date, ii);
    
    suptitle(h,filename,1);
    
    % saving as png
    savepngFast(h,rd.savingpath,pngfilename);
    close(h);
    % update message
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', rd.clusters(ii,1),rd.clusters(ii,2));
    fprintf(msg);
    
end

disp(seconds2human (toc(totaltime)));
sound(struct2array(load('chirp.mat','y')))


end
