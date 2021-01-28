

function plotDSdata(ds, clus, varargin)
%
%%% plotDS %%%
%
%
% This function plots the polar plot for showing the preferred direction of
% the response.
%
%================================Inputs====================================
%
%   nbins : number of bins for the histograms
%   nangap : gap between the cycles of raster plots
%   colords : DS plot color.
%   colordslines : DS plot line colors.
%   position : figure size [x y width height].
%   savefigure : flag for saving the figure.
%   savedata : flag for saving the data for each individual cell.
%   figuresavepath : figure saving path.
%   figurname : main title that is used for figure joined with rgcname function.
%   visible : figure visibility.
%   title : starting title to separate moving bar from grating stimuli.
%   savealltogether : save everthing in one file in AllData format.
%
%================================Output====================================
%
%   DS plot : the main output is DS plot.
%
% written by Mohammad, 14.02.2018.
% major update to faster version with new plotting for CMOS analysis on
% 11.05.2019.
% minor update for the opto project on 13.01.2021.

% first parse the inputs
p = inputParser();
p.addParameter('nbins', 50, @(x) isnumeric(x));
p.addParameter('nangap', 2, @(x) isnumeric(x));
p.addParameter('colords', rgb('azure'), @(x) isnumeric(x));
p.addParameter('colordslines', rgb('deepskyblue'), @(x) isnumeric(x));
p.addParameter('position', [20 500-(170*size(ds,2)) 1800 350*size(ds,2)], @(x) isnumeric(x));
p.addParameter('savefigure', false, @(x) islogical(x) );
p.addParameter('savedata', false, @(x) islogical(x) );
p.addParameter('figuresavepath', [], @(x) ischar(x));
p.addParameter('datasavepath', [], @(x) ischar(x));
p.addParameter('figurname', 'Direction Grating Sequence', @(x) ischar(x));
p.addParameter('visible', 'off', @(x) ischar(x));
p.addParameter('title', [], @(x) ischar(x));
p.addParameter('savealltogether', false, @(x) islogical(x));
% p.addParameter('oldversion', false, @(x) islogical(x));
p.addParameter('sortinfo', []);
p.parse(varargin{:});
pltpara = p.Results;

if pltpara.savefigure && isempty(pltpara.figuresavepath)
    error('Aint nobody got a path for saving figure, put figuresavepath in the input');
end

if pltpara.savedata && isempty(pltpara.datasavepath)
    error('Aint nobody got a path for saving data, put figuresavepath in the input');
end

% first calculate everything then plot
nang = [ds(1,:).para];
nang = [nang.nangles];
cira = circAvgallDScells(ds, nang, size(nang,2));
hstdist = getdistfit(ds, pltpara.nbins);
rasters = preprasters(ds,pltpara.nangap);

numcells = size(ds,1);
numstim = size(ds,2);
col = lines(ds(1,1).para.nangles);
angs = linspace(0,360,(ds(1,1).para.nangles)+1); angs = angs(1:end-1);
dsiall = reshape([ds.dsi],[],size(ds,2));
osiall = reshape([ds.osi],[],size(ds,2));
dsipvalall = reshape([ds.dsi_pval],[],size(ds,2));
osipvalall = reshape([ds.osi_pval],[],size(ds,2));
msg=[];
clearvars p nang;
warning('off','MATLAB:prnRenderer:opengl');
for ii = 1:numcells
    
    if ii == 1
        tic;
        h = figure('pos',pltpara.position,'color','w','visible',pltpara.visible);
        p = cell(numstim,1);      [ax1,ax2] = deal(zeros(1,numstim));
        for jj = 1: numstim
            dsseqdata = ds(ii,jj);
            para = dsseqdata.para;
            dsp = subplot_tight(numstim,6,1+(jj-1)*6,0.035);
            p{jj}.dsgr = plot(squeeze(cira.grx(ii,jj,:)),squeeze(cira.gry(ii,jj,:)),'color',rgb('gray'),'LineWidth',0.25);
            hold on;
            p{jj}.dsp = plot(squeeze(cira.x(ii,jj,:)),squeeze(cira.y(ii,jj,:)),'color',pltpara.colordslines,'LineWidth',2 );
            p{jj}.cirv = plot([0 cira.circAvg(ii,jj,1)],[0 cira.circAvg(ii,jj,2)],'color',pltpara.colordslines,'LineWidth',2);
            p{jj}.cirvspot = plot(cira.circAvg(ii,jj,1),cira.circAvg(ii,jj,2),'o','color',pltpara.colordslines,...
                'LineWidth',3,'MarkerFaceColor',pltpara.colordslines,'MarkerSize',6);
            p{jj}.dstxt1 = text( cira.txtvals(ii,jj,1),0, {sprintf(' %g', round(cira.txtvals(ii,jj,1)));' Hz'},'verticalalignment', ...
                'middle','horizontalalignment', 'left','fontsize',7);
            p{jj}.dstxt2 = text( cira.txtvals(ii,jj,2),0, {sprintf(' %g', round(cira.txtvals(ii,jj,2)));' Hz'},'verticalalignment', ...
                'middle','horizontalalignment', 'left','fontsize',7);
            axis equal tight off;
            dsppos = get(dsp,'pos'); set(dsp,'pos',[dsppos(1)-0.02, dsppos(2),dsppos(3), dsppos(4)]);
            
            % rasters
            subplot_tight(size(ds,2),6,2+(jj-1)*6,0.035)
            
            yAx = (pltpara.nangap+0.5:(para.cycles+pltpara.nangap):(para.cycles+pltpara.nangap)*para.nangles);
            
            if para.regeneration==0
                bkgdur = 0;
            else
                bkgdur = para.regeneration/4/para.fps; % the last quarter of the gray frames are considered
                if bkgdur <= 0.1, bkgdur = 0.25; end  % cannot be shorter that  250 ms
            end
            pdur = bkgdur;
            firsttrialdur =  para.period/para.fps;
            rectangle('pos',[0, 0, pdur max(yAx)+pltpara.nangap],'facecolor',rgb('gray 90'),'edgecolor','none'); hold on;
            rectangle('pos',[pdur, 0, firsttrialdur max(yAx)+pltpara.nangap],'facecolor',rgb('lightpink'),'edgecolor','none');
            %maxX = zeros(1,para.nangles);
            for kk = 1:para.nangles
                try
                    p{jj}.pras(kk) = line(squeeze(rasters(ii,jj,kk).x),squeeze(rasters(ii,jj,kk).y),'color',col(kk,:),'LineWidth',0.35);
                catch
                    p{jj}.pras(kk) = line(NaN,NaN, 'color',col(kk,:),'LineWidth',0.35);
                end
                try
                    p{jj}.pfrras(kk) = line(squeeze(rasters(ii,jj,kk).pfrx),squeeze(rasters(ii,jj,kk).pfry),'color',0.25*[1 1 1],'LineWidth',0.35);
                catch
                    p{jj}.pfrras(kk) = line(NaN,NaN, 'color',0.25*[1 1 1],'LineWidth',0.35);
                end
                try
                    p{jj}.firstras(kk) = line(squeeze(rasters(ii,jj,kk).trial1x),squeeze(rasters(ii,jj,kk).trial1y),'color','k','LineWidth',0.35);
                catch
                    p{jj}.firstras(kk) = line(NaN,NaN,'color','k','LineWidth',0.35);
                end
            end
            xAx = pdur + (para.duration/para.fps) ;     if xAx < 0.1 || isnan(xAx), xAx = 5; end
            axis([0 xAx, 0 max(yAx)+pltpara.nangap]);   axis square;
            ylabel('Angles'); xlabel('Time (sec)');    if jj==1, title('Rasters','fontsize',8); end
            xticks([0 pdur pdur+firsttrialdur:firsttrialdur:xAx]);
            xticklabels(round([0 pdur pdur+firsttrialdur:firsttrialdur:xAx],1)-round(pdur,2));
            set(gca,'ytick',yAx,'yticklabel',angs,'fontsize',7,'ticklength',[0.015, 0.015]);
            if numel(para.period) > 3, pbaspect([2,1,1]); else, pbaspect([1,1,1]); end
            
            % dsi p-values compared to the shuffeled analysis
            ax1(jj) = subplot_tight(size(ds,2),6,3+(jj-1)*6,0.035);
            p{jj}.hdsb = bar(squeeze (hstdist.dshistx(ii,jj,:)),squeeze(hstdist.dshist(ii,jj,:)),'facecolor',rgb('sandybrown'),'edgecolor','none','BarWidth',1);
            p{jj}.hdsl = line(squeeze(hstdist.dsx(ii,jj,:)),squeeze(hstdist.dsnormfit(ii,jj,:)),'color','k','linewidth',2);
            hold on;        box off; axis square;
            p{jj}.hdss = stem(dsseqdata.dsi,10,'o','markersize',7,'markerfacecolor','r','color',rgb('crimson'),'linewidth',2);
            p{jj}.hdsleg = legend(p{jj}.hdss,['DSI: ',num2str(round(dsseqdata.dsi,2)),newline,'p-value: ',num2str(round(dsseqdata.dsi_pval,3))]);
            %lb(2).Children.Children(2).LineStyle = 'none';
            legend boxoff;
            yAx = ceil(max(squeeze(hstdist.dsnormfit(ii,jj,:)))/2)*2;               if yAx <1 || isnan(yAx), yAx = 5; end
            xAx = max([max(squeeze(hstdist.dsx(ii,jj,:))),dsseqdata.dsi])+0.05;     if xAx <=0|| isnan(xAx), xAx = 1; end
            axis([-0.015 xAx  0 yAx]);
            ylabel('Shffeled bins'); xlabel('DSI');     set(gca,'ytick',0:yAx/2:yAx,'fontsize',7,'ticklength',[0.015, 0.015]);
            if jj==1, title('DSI vs shuffeled','fontsize',8); end
            if numel(para.period) > 3, pbaspect([2,1,1]); else, pbaspect([1,1,1]); end
            
            
            % osi p-values compared to the shuffeled analysis
            ax2(jj) = subplot_tight(size(ds,2),6,4+(jj-1)*6,0.035);
            p{jj}.hosb = bar(squeeze (hstdist.oshistx(ii,jj,:)),squeeze(hstdist.oshist(ii,jj,:)),'facecolor',rgb('skyblue'),'edgecolor','none','BarWidth',1);
            p{jj}.hosl = line(squeeze(hstdist.osx(ii,jj,:)),squeeze(hstdist.osnormfit(ii,jj,:)),'color','k','linewidth',2);
            hold on;        box off; axis square;
            p{jj}.hoss = stem(dsseqdata.osi,10,'o','markersize',7,'markerfacecolor',rgb('violet'),'color',rgb('purple 1'),'linewidth',2);
            p{jj}.hosleg = legend(p{jj}.hoss,['OSI: ',num2str(round(dsseqdata.osi,2)),newline,'p-value: ',num2str(round(dsseqdata.osi_pval,3))]);
            %lb(2).Children.Children(2).LineStyle = 'none';
            legend boxoff;
            yAx = ceil(max(squeeze(hstdist.osnormfit(ii,jj,:)))/2)*2;               if yAx <1 || isnan(yAx), yAx = 5; end
            xAx = max([max(squeeze(hstdist.osx(ii,jj,:))),dsseqdata.osi])+0.05;     if xAx <=0 || isnan(xAx), xAx = 1; end
            axis([-0.015 xAx 0 yAx]);
            ylabel('Shffeled bins'); xlabel('OSI');     set(gca,'ytick',0:yAx/2:yAx,'fontsize',7,'ticklength',[0.015, 0.015]);
            if jj==1, title('OSI vs shuffeled','fontsize',8); end
            if numel(para.period) > 3, pbaspect([2,1,1]); else, pbaspect([1,1,1]); end
            
            
            % dsi compared to the osi
            subplot_tight(size(ds,2),6,5+(jj-1)*6,0.035)
            plot(dsiall(:,jj),osiall(:,jj),'o','markerfacecolor',rgb('silver'),'color',rgb('gray'),'markersize',3);
            hold on;        plot([0 1],[0 1],'k--'); axis equal square;
            yAxdo = round(max([max(dsiall(:,jj)),max(osiall(:,jj))]),3);          if yAxdo <0.1 || isnan(yAxdo), yAxdo = 0.5; end
            axis([0 yAxdo+0.01 0 yAxdo+0.01]);      box off;
            p{jj}.pdsos = plot(dsseqdata.dsi,dsseqdata.osi,'o','markerfacecolor',rgb('red'),'color',rgb('crimson'),'markersize',8);
            ylabel('OSI'); xlabel('DSI');     set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1,'fontsize',7,'ticklength',[0.02, 0.02]);
            if jj==1, title('DSI vs OSI','fontsize',8); end
            if numel(para.period) > 3, pbaspect([2,1,1]); else, pbaspect([1,1,1]); end
            
            
            % p-values for the dsi and osi
            subplot_tight(size(ds,2),6,6+(jj-1)*6,0.035)
            plot(dsipvalall(:,jj),dsiall(:,jj),'o','markerfacecolor','none','color',rgb('salmon'),'markersize',4); hold on;
            plot(osipvalall(:,jj),osiall(:,jj),'o','markerfacecolor','none','color',rgb('skyblue'),'markersize',4);
            axis([5e-4 1.5 -0.01 yAxdo+0.01]);    axis square; set(gca,'xscale','log');
            p{jj}.pdsi = plot(dsseqdata.dsi_pval,dsseqdata.dsi,'o','markerfacecolor',rgb('red'),'color',rgb('black'),'markersize',8,'linewidth',2);
            p{jj}.posi = plot(dsseqdata.osi_pval,dsseqdata.osi,'o','markerfacecolor',rgb('royalblue'),'color',rgb('black'),'markersize',8,'linewidth',2);
            if jj==1, title('p-values for DSI & OSI','fontsize',8); end
            ylabel('DSI & OSI');    xlabel('p-values');         box off;
            set(gca,'xscale','log','xtick',[1e-3,1e-2,1e-1,1],'ytick',0:0.2:1,'fontsize',7,'ticklength',[0.02, 0.02]);
            if numel(para.period) > 3, pbaspect([2,1,1]); else, pbaspect([1,1,1]); end
            
        end
        
        fprintf('First plot time: ''%g'' sec\n',round(toc,3));
    else
        
        for jj = 1: numstim
            dsseqdata = ds(ii,jj);
            % ds plot
            p{jj}.dsgr.XData = squeeze(cira.grx(ii,jj,:));      p{jj}.dsgr.YData = squeeze(cira.gry(ii,jj,:));
            p{jj}.dsp.XData = squeeze(cira.x(ii,jj,:));         p{jj}.dsp.YData = squeeze(cira.y(ii,jj,:));
            p{jj}.cirv.XData = [0 cira.circAvg(ii,jj,1)];       p{jj}.cirv.YData = [0 cira.circAvg(ii,jj,2)];
            p{jj}.cirvspot.XData = cira.circAvg(ii,jj,1);       p{jj}.cirvspot.YData = cira.circAvg(ii,jj,2);
            p{jj}.dstxt1.Position = [cira.txtvals(ii,jj,1),0,0];
            p{jj}.dstxt1.String = {sprintf(' %g', round(cira.txtvals(ii,jj,1)));' Hz'};
            p{jj}.dstxt2.Position = [cira.txtvals(ii,jj,2),0,0];
            p{jj}.dstxt2.String = {sprintf(' %g', round(cira.txtvals(ii,jj,2)));' Hz'};
            
            %rasters plots
            for kk = 1:para.nangles
                p{jj}.pras(kk).XData = squeeze(rasters(ii,jj,kk).x);    p{jj}.pras(kk).YData = squeeze(rasters(ii,jj,kk).y);
                p{jj}.pfrras(kk).XData = squeeze(rasters(ii,jj,kk).pfrx);   p{jj}.pfrras(kk).YData = squeeze(rasters(ii,jj,kk).pfry);
                p{jj}.firstras(kk).XData = squeeze(rasters(ii,jj,kk).trial1x);
                p{jj}.firstras(kk).YData = squeeze(rasters(ii,jj,kk).trial1y);
            end
            % ds hist
            p{jj}.hdsb.XData = squeeze (hstdist.dshistx(ii,jj,:));   p{jj}.hdsb.YData = squeeze(hstdist.dshist(ii,jj,:));
            p{jj}.hdsl.XData = hstdist.dsx(ii,jj,:);                 p{jj}.hdsl.YData = squeeze(hstdist.dsnormfit(ii,jj,:));
            p{jj}.hdss.XData = dsseqdata.dsi;
            p{jj}.hdsleg.String = ['dsi: ',num2str(round(dsseqdata.dsi,2)),newline,'p-value: ',num2str(round(dsseqdata.dsi_pval,3))];
            yAx = ceil(max(squeeze(hstdist.dsnormfit(ii,jj,:)))/2)*2;               if yAx <=0 || isnan(yAx), yAx = 5; end
            xAx = max([max(squeeze(hstdist.dsx(ii,jj,:))),dsseqdata.dsi])+0.05;     if xAx <=0 || isnan(xAx), xAx = 1; end
            axis(ax1(jj), [-0.015 xAx  0 yAx]);            set(ax1(jj),'ytick',0:yAx/2:yAx);
            % os hist
            p{jj}.hosb.XData = squeeze (hstdist.oshistx(ii,jj,:));   p{jj}.hosb.YData = squeeze(hstdist.oshist(ii,jj,:));
            p{jj}.hosl.XData = hstdist.osx(ii,jj,:);                 p{jj}.hosl.YData = squeeze(hstdist.osnormfit(ii,jj,:));
            p{jj}.hoss.XData = dsseqdata.osi;
            p{jj}.hosleg.String = ['dsi: ',num2str(round(dsseqdata.osi,2)),newline,'p-value: ',num2str(round(dsseqdata.osi_pval,3))];
            yAx = ceil(max(squeeze(hstdist.osnormfit(ii,jj,:)))/2)*2;               if yAx <=0 || isnan(yAx), yAx = 5; end
            xAx = max([max(squeeze(hstdist.osx(ii,jj,:))),dsseqdata.osi])+0.05;     if xAx <=0 || isnan(xAx), xAx = 1; end
            axis(ax2(jj),[-0.015 xAx 0 yAx]);            set(ax2(jj),'ytick',0:yAx/2:yAx);
            
            % dsos
            p{jj}.pdsos.XData = dsseqdata.dsi;          p{jj}.pdsos.YData = dsseqdata.osi;
            % ds os p-values
            p{jj}.pdsi.XData = dsseqdata.dsi_pval;      p{jj}.pdsi.YData = dsseqdata.dsi;
            p{jj}.posi.XData = dsseqdata.osi_pval;      p{jj}.posi.YData = dsseqdata.osi;
        end
    end
    
    if isempty(pltpara.title) && ~isempty(pltpara.datasavepath)
        
        if isstruct(pltpara.sortinfo), chinfo =  pltpara.sortinfo(ii); else, chinfo = clus(ii,:); end
        [filename,pngfilename] = rgcname(pltpara.figurname, chinfo, para.date, ii);
    else
        filename = pltpara.title;
        if isempty(filename), filename = 'Direction Grating Sequence'; end
    end
    
    if ii==1
        sm = suptitle(h,filename,2);
        smorigpos = sm.Position;
    else
        sm.String = filename;
        sm.Position = smorigpos;
    end
    
    if pltpara.savefigure
        savepngFast(h,pltpara.figuresavepath,pngfilename); % a bit lower quality for faster running
    end
    
    if pltpara.savedata && ~pltpara.savealltogether
        switch lower(pltpara.figurname)
            case 'direction grating sequence'
                dsseqdata = ds(ii,:);
                save([pltpara.datasavepath,pngfilename],'dsseqdata');
            case 'moving bar'
                mbdata = dsseqdata;
                save([pltpara.datasavepath,pngfilename],'-struct','mbdata');
            case 'direction grating'
                dsdata = ds(ii,:);
                save([pltpara.datasavepath,pngfilename],'dsdata');
        end
    end
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', clus(ii,1),clus(ii,2));
    fprintf(msg);
    %disp(['Analysis for Cell ',num2str(clus(ii,1)),', Cluster ',num2str(clus(ii,2)),' is... wait for it...done!!!']);
end
if pltpara.savealltogether
    [~,allsavename] = fileparts(pltpara.figuresavepath(1:end-1));
    allsavename = strrep([allsavename,'_for_experiment_on_',para.date],'Analysis','analysis');
    save([pltpara.datasavepath,allsavename],'ds','-v7.3');
end
close(h);
warning('on','MATLAB:prnRenderer:opengl');
end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function hs = getdistfit(ds, nbins)

dsdist = cell2mat({ds.dsi_dist}');
osdist = cell2mat({ds.osi_dist}');
tic;
[dsh, dshx, osh, oshx] = deal(nan(size(dsdist,1),nbins));
[dsnormfit, dsx, osnormfit, osx] = deal(nan(size(dsdist,1),100));
for ii = 1:size(dsdist,1)
    
    if  ~all(isnan(dsdist(ii,:)))
        % these are duplicate from histfit function
        pd = fitdist(dsdist(ii,:)','normal');
        % Find range for plotting
        q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
        dsx(ii,:) = linspace(q(1),q(2),100);
        
        % Do histogram calculations
        [dsh(ii,:),binedges] = histcounts(dsdist(ii,:),nbins);
        dshx(ii,:) = binedges(1:end-1)+diff(binedges)/2;
        % Normalize the density to match the total area of the histogram
        binwidth = binedges(2)-binedges(1); % Finds the width of each bin
        area = size(dsdist,2) * binwidth;
        dsnormfit(ii,:) = area * pdf(pd, dsx(ii,:));
    end
    %--------------OS cells---------------------------------------------
    if  ~all(isnan(osdist(ii,:)))
        pd = fitdist(osdist(ii,:)','normal');
        % Find range for plotting
        q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
        osx(ii,:) = linspace(q(1),q(2),100);
        
        % Do histogram calculations
        [osh(ii,:),binedges] = histcounts(dsdist(ii,:),nbins);
        oshx(ii,:) = binedges(1:end-1)+diff(binedges)/2;
        % Normalize the density to match the total area of the histogram
        binwidth = binedges(2)-binedges(1); % Finds the width of each bin
        area = size(osdist,2) * binwidth;
        osnormfit(ii,:) = area * pdf(pd, dsx(ii,:));
    end
end

hs.dshist = reshape(dsh,size(ds,1),size(ds,2),[]);
hs.dshistx = reshape(dshx,size(ds,1),size(ds,2),[]);
hs.dsx = reshape(dsx,size(ds,1),size(ds,2),[]);
hs.dsnormfit = reshape(dsnormfit,size(ds,1),size(ds,2),[]);

hs.oshist = reshape(osh,size(ds,1),size(ds,2),[]);
hs.oshistx = reshape(oshx,size(ds,1),size(ds,2),[]);
hs.osx = reshape(osx,size(ds,1),size(ds,2),[]);
hs.osnormfit = reshape(osnormfit,size(ds,1),size(ds,2),[]);

fprintf('histogram calculation time: ''%g'' sec\n',round(toc,3));

end

%--------------------------------------------------------------------------------------------------%

function ras = preprasters(ds,nangap)

tic;
spkheight = 0.4;
ras = cell([size(ds),8]);

for ii = 1:size(ds,1)
    for jj = 1: size(ds,2)
        
        if ds(ii,jj).para.regeneration==0
            bkgdur = 0;
        else
            bkgdur = ds(ii,jj).para.regeneration/4/ds(ii,jj).para.fps; % the last quarter of the gray frames are considered
            if bkgdur <= 0.1, bkgdur = 0.25; end  % cannot be shorter that  250 ms
        end
        pdur = bkgdur;
        firsttrialdur =  ds(ii,jj).para.period/ds(ii,jj).para.fps;
        rasloc = (0:ds(ii,jj).para.nangles-1)* (ds(ii,jj).para.cycles+nangap);
        
        for kk = 1:ds(ii,jj).para.nangles
            
            rasresp = firsttrialdur + pdur + ds(ii,jj).rasters{kk};
            spkloc = repmat(1:size(rasresp,1),size(rasresp,2),1)';
            spkoneline = [rasresp(:),rasresp(:),nan(numel(rasresp(:)),1)]';
            tboneline = [spkloc(:)-spkheight,spkloc(:)+spkheight, nan(numel(spkloc(:)),1)]' + rasloc(kk) ;
            
            ras{ii,jj,kk}.x = spkoneline(:);
            ras{ii,jj,kk}.y = tboneline(:);
            
            pfrras = ds(ii,jj).pfrrasters{kk};
            spkloc = repmat(1:size(pfrras,1),size(pfrras,2),1)';
            spkoneline = [pfrras(:),pfrras(:),nan(numel(pfrras(:)),1)]';
            tboneline = [spkloc(:)-spkheight,spkloc(:)+spkheight, nan(numel(spkloc(:)),1)]' + rasloc(kk) ;
            
            ras{ii,jj,kk}.pfrx = spkoneline(:);
            ras{ii,jj,kk}.pfry = tboneline(:);
            
            frtrras = pdur+ds(ii,jj).firsttrialrasters{kk};
            spkloc = repmat(1:size(frtrras,1),size(frtrras,2),1)';
            spkoneline = [frtrras(:),frtrras(:),nan(numel(frtrras(:)),1)]';
            tboneline = [spkloc(:)-spkheight,spkloc(:)+spkheight, nan(numel(spkloc(:)),1)]' + rasloc(kk);
            
            ras{ii,jj,kk}.trial1x = spkoneline(:);
            ras{ii,jj,kk}.trial1y = tboneline(:);
            
            clearvars nanras pfrnanras frtrnanras spkloc spkoneline tboneline;
        end
        
    end
end
ras = cell2mat(ras);
fprintf('rasters calculation time: ''%g'' sec\n',round(toc,3));
end
