

function [] = plot_Marmoset_ReversingGratingWVSP(datapath, varargin)
%ANALYZE_GRATINGS Analyzes responses to OnOffGratings stimuli
%Also plots analyses for each cell

dp = dir([datapath,filesep,'Data Analysis/','*ReversingGratingWVSP_Analysis*']);
dp = [dp.folder,filesep,dp.name];
%load data
%----------------------------------------------------------------------
rdata = dir([dp,filesep,'rgwvsp_data/*All reversinggratingwvsp_analysis data*.mat']);
rdata = load([rdata.folder,filesep,rdata.name]);

clus = struct2array(load([datapath,'/',findFileinFolder(datapath,'CellsList','mat')]));
% clus = clus(clus(:,3)<=2,:);
stimPara=rdata.para;
%----------------------------------------------------------------------
Ncells   = size(clus,1);
maxPhases = max(stimPara.nphases);
Nwidths = numel(stimPara.stripewidth);
Nreversals = size(rdata.allRasters,3);
Ntrials = size(rdata.allRasters,5);
sphases = stimPara.nphases;
%xdata = [0 stimPara.stripewidth];
% if xdata(end)>=stimPara.screen(1)
%     ff_flag=1; xdata(end)=[];
% else ff_flag=0;
% end
swidths = stimPara.stripewidth * stimPara.monitorpixel; % experiment.projector.pixelsize * 1e6;
cols = flipud(colormap(cbrewer('seq','YlGnBu',Ntrials+2))); % linspace(0,0.4,Ntrials)'.*[1 1 1];
xfreq = 1e3./(2*swidths);
close all;
%----------------------------------------------------------------------
tic;
h = figure('visible','off','Position', [150 10 200*maxPhases 120*Nwidths],'Color','w');
h.Units = 'centimeters';
%----------------------------------------------------------------------
ptimes = rdata.allTimes;
ptimes = repmat(ptimes,[2 1])+[-1;1]*mean(diff(ptimes))/2;
ptimes = [ptimes(1); ptimes(:); ptimes(end)];
%----------------------------------------------------------------------
p = panel();
p.pack('h',{0.75 0.25}); p(1).pack('v',Nwidths);

p(2).pack('v',2);

for iwidth = 1:Nwidths
    p(1,iwidth).pack('h',maxPhases);
    for iphase = 1:maxPhases
        p(1,iwidth,iphase).pack('v',{0.6 0.4});
    end
end

p.de.margin = 2; p.fontsize = 10;
p(1).de.marginleft = 5;
p(2).marginleft = 20; p(2).de.margintop= 30;
p.margin = [12 12 8 12];
%----------------------------------------------------------------------
for iwidth = 1:Nwidths
    for iphase = 1:maxPhases
        p(1,iwidth,iphase,1).select(); ax=gca; ax.Visible='off';
        ax.YDir = 'reverse';
        ylim(0.5+ [0 Ntrials*Nreversals]); xlim([1 rdata.trialdur]);
        
        p(1,iwidth,iphase,2).select(); ax=gca; ax.Visible='off';
        ylim([0 1]); xlim([min(ptimes) max(ptimes)])
        p(1,iwidth).ylabel(sprintf('%1.1f µm', swidths(iwidth)));
    end
end
%----------------------------------------------------------------------
p(2,1).select(); axis square; xlim([0 max(swidths)])
xlabel('Bar width (µm)');       ylabel('Firing rate (Hz)');
p(2,2).select(); axis square; xlim([min(xfreq) max(xfreq)]);
xlabel('Spatial frequency (cyc/mm)');       ylabel('Amplitude');
%----------------------------------------------------------------------
spkheight = 0.4;
msg=[];
for icell = 1:Ncells    
    % ----------------------------------------------------------------------
    cellRasters = squeeze(rdata.allRasters(icell,:,:,:,:));
    cellRates   = squeeze(rdata.allRates(icell,:,:));
    maxr = max([max(cellRates,[],'all') 20]);
    %----------------------------------------------------------------------
    for iw = 1:Nwidths
        
        winds = sum(sphases(1:iw-1))+1:sum(sphases(1:iw-1))+sphases(iw);
        for ip = 1:numel(winds)
            plotRaster = squeeze(logical(full(cellRasters(:,:,winds(ip),:))));
            
            p(1,iw,ip,1).select(); cla;
            for it = 1:Ntrials                 
                [rx,ry] = find(plotRaster(:,:,it));
                %rx = rdata.rasterTimes(rx)*1e3;
                spkxline = [rx(:),rx(:),nan(numel(rx(:)),1)]';
                spkyline = [ry(:)-spkheight,ry(:)+spkheight, nan(numel(ry(:)),1)]' ;
                line(spkxline(:),spkyline(:)+Nreversals*(it-1),'color',cols(it,:),'linewidth',0.5);
                %rasterPlotter(plotRaster(:,:,it)', 0.4, cols(it,:), 0.5, Nreversals*(it-1));
            end            
            p(1,iw,ip,2).select(); cla;
            frates = cellRates(:,winds(ip))/maxr;
            patch(ptimes,[0;kron(frates,ones(2,1));0],rgb('deepskyblue'),'EdgeColor',rgb('navy'),'linewidth',0.025);
        end
    end
    %----------------------------------------------------------------------
    p(2,1).select();
    ydata = rdata.maxRates(icell,:);
    yyfit = logistic4(rdata.allParams(icell,:),[0 swidths]);
    ylim([0 maxr]);
    semilogx([0 swidths], ydata, '-ok', [0 swidths], yyfit, '-or','LineWidth',2);
    title(['Quality Rsq: ',num2str(round(rdata.qualityRsq(icell),2))]);
    %----------------------------------------------------------------------
    p(2,2).select(); cla;
    f1comp = rdata.allF1(icell,:);
    f2comp = rdata.allF2mean(icell,:);
    semilogx(xfreq, f1comp, '-ok',xfreq, f2comp, '-or', 'LineWidth',2);
    legend('F1','F2'); legend boxoff; 
    title(['Nonlin index: ',num2str(round(rdata.nlinIdx(icell),2)),...
        ' , Classic nonlin index: ',num2str(round(rdata.classicNlinIdx(icell),2))]);
    
    if isfield(rdata,'sortinfo')
        [filename,filenamesavepng] = generateRGCname('Reversing Grating WVSP', rdata.sortinfo(icell),rdata.savingpath);
    else
        filename = generateRGCname('Reversing Grating WVSP', clus(icell,:),rdata.savingpath);
    end    
    % this is to have shorter name for png saving
    % this is to have shorter name for png saving
    if exist('filenamesavepng','var')
        pngfilename = [num2str(icell,'%02d-'),extractAfter(filenamesavepng,[rdata.savingpath,'\'])];
    else
        pngfilename = filename;
    end
    %filename = generateRGCname('Reversing Grating WVSP',clus(icell,:),rdata.savingpath);
    p.title(filename);
    savepngFast(h,rdata.savingpath, pngfilename);
    
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', icell,Ncells,toc);
    fprintf(msg);
end
close(h);
%----------------------------------------------------------------------

end


