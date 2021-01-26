

function rgdata = plot_ReversingGratingWVSP(datapath, varargin)
%
%%% plot_ReversingGratingWVSP %%%
%
%
% This function plot the analyzed data for the reversing grating with varying
% spatial period (WVSP) stimulus.
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder. It can also be the analyzed data
%              struture from the Analyze_ReversingGrating function.
%
%================================Output====================================
%
%   rgdata : output is the structure of the analyzed data.
%   Plot : Rasters and PSTHs for all the grating width and spatial periods.
%          Also the F1 and F2 plot from the responses.
%
% written by Mohammad based on similar function from Dimos.
% moved to compeletly new version on for the opto project 12.01.2021.

%load data
%----------------------------------------------------------------------
if isstruct(datapath)
    rgdata = datapath;
else
    rgdp = dir([datapath,filesep,'Data Analysis',filesep,'*ReversingGratingWVSP_Analysis*']);
    
    if size(rgdp,1) > 1
        selectmsg = 'Too many reversing grating stimuli on the dance floor!! select what you want to analyze!! ==> ';
        for ii=1:   size(rgdp,1), disp(rgdp(ii).name);        end
        expnum = input(selectmsg);
        rgdp = rgdp(expnum);
    end
    
    rgdatpath = dir([[rgdp.folder,filesep,rgdp.name,filesep],'*Reversing_grating_wvsp_analysis*']);
    if size(rgdatpath,1) > 1
        selectmsg = 'Too many reversing grating data on the dance floor!! select what you want to analyze!! ==> ';
        for ii=1:   size(rgdatpath,1), disp(rgdatpath(ii).name);        end
        expnum = input(selectmsg);
        rgdatpath = rgdatpath(expnum);
    end
    rgdata = load([rgdatpath.folder,filesep,rgdatpath.name]);
end

stimPara    =   rgdata.para;
Ncells      =   size(rgdata.clusters,1);
maxPhases   =   max(stimPara.Nphases);
Nwidths     =   numel(stimPara.stripewidths);
Nreversals  =   size(rgdata.allRasters,3);
Ntrials     =   size(rgdata.allRasters,5);
sphases     =   stimPara.Nphases;
swidths     =   stimPara.stripewidths * stimPara.pixelsize; % experiment.projector.pixelsize * 1e6;
cols        =   lines(Nwidths);
xfreq       =   1e3./(2*swidths);
%----------------------------------------------------------------------
tic;
h = figure('visible','off','Position', [150 10 200*maxPhases 120*Nwidths],'Color','w');
h.Units = 'centimeters';
%----------------------------------------------------------------------
ptimes = rgdata.allTimes;
ptimes = repmat(ptimes,[2 1])+[-1;1]*mean(diff(ptimes))/2;
ptimes = [ptimes(1); ptimes(:); ptimes(end)];
%----------------------------------------------------------------------
p = panel();
p.pack('h',{0.8 0.2}); p(1).pack('v',Nwidths);

p(2).pack('v',2);

for iwidth = 1:Nwidths
    p(1,iwidth).pack('h',maxPhases);
    for iphase = 1:maxPhases
        p(1,iwidth,iphase).pack('v',{0.6 0.4});
    end
end

p.de.margin = 2; p.fontsize = 11;
p(1).de.marginleft = 6;
p(1).de.margintop = 2;
p(2).marginleft = 20; p(2).de.margintop= 2;
p.margin = [15 2 10 10];
%----------------------------------------------------------------------
for iwidth = 1:Nwidths
    for iphase = 1:maxPhases
        p(1,iwidth,iphase,1).select(); ax=gca; ax.Visible='off';
        ax.YDir = 'reverse';
        ylim(0.5+ [0 Ntrials*Nreversals]); xlim([1 rgdata.trialdur]);
        
        p(1,iwidth,iphase,2).select(); ax=gca; ax.Visible='off';
        ylim([0 1]); xlim([min(ptimes) max(ptimes)])
        p(1,iwidth).ylabel(sprintf('%1.1f µm', swidths(iwidth)));
    end
end
%----------------------------------------------------------------------
p(2,1).select(); axis square; xlim([-1 max(swidths)])
xlabel('Bar width (µm)');       ylabel('Firing rate (Hz)');
xticks(logspace(-4,10,15));
p(2,2).select(); axis square; xlim([min([-1,xfreq]) max(xfreq)]);
xlabel('Spatial frequency (cyc/mm)');       ylabel('Amplitude');
xticks(logspace(-4,10,15));
%----------------------------------------------------------------------
spkheight = 0.4;
msg=[];

for icell = 1:Ncells
    % ----------------------------------------------------------------------
    cellRasters = squeeze(rgdata.allRasters(icell,:,:,:,:));
    cellRates   = squeeze(rgdata.allRates(icell,:,:));
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
                line(spkxline(:),spkyline(:)+Nreversals*(it-1),'Color',cols(iw,:),'LineWidth',0.5);
                %rasterPlotter(plotRaster(:,:,it)', 0.4, cols(it,:), 0.5, Nreversals*(it-1));
            end
            p(1,iw,ip,2).select(); cla;
            frates = cellRates(:,winds(ip))/maxr;
            patch(ptimes,[0;kron(frates,ones(2,1));0],cols(iw,:),'EdgeColor','k','LineWidth',0.025);
        end
    end
    %----------------------------------------------------------------------
    p(2,1).select(); cla;
    ydata = rgdata.maxRates(icell,:);
    yyfit = logistic4(rgdata.allParams(icell,:),[0 swidths]);
    ylim([0 maxr+2]);
    semilogx([0 swidths], ydata, '-o','Color',cols(1,:),'MarkerFaceColor',cols(6,:), 'LineWidth',2);
    hold on;
    semilogx([0 swidths], yyfit, '-or','Color',cols(7,:),'MarkerFaceColor',cols(2,:),'LineWidth',2);
    title(['Quality Rsq: ',num2str(round(rgdata.qualityRsq(icell),2))]);
    %----------------------------------------------------------------------
    p(2,2).select(); cla;
    f1comp = rgdata.allF1(icell,:);
    f2comp = rgdata.allF2mean(icell,:);
    semilogx(xfreq, f1comp, '-o','Color',cols(1,:),'MarkerFaceColor',cols(6,:), 'LineWidth',2);
    hold on;
    semilogx(xfreq, f2comp, '-or','Color',cols(7,:),'MarkerFaceColor',cols(2,:),'LineWidth',2);
    %semilogx(xfreq, f1comp, '-ok',xfreq, f2comp, '-or', 'LineWidth',2);
    legend('F1','F2'); legend boxoff;
    title(['Nonlin index: ',num2str(round(rgdata.nlinIdx(icell),2)),...
        ' , Classic nonlin index: ',num2str(round(rgdata.classicNlinIdx(icell),2))]);
    
    if isfield(rgdata,'sortinfo')
        [filename,filenamesavepng] = rgcname('Reversing Grating WVSP', rgdata.sortinfo(icell),stimPara.date);
        pngfilename = [num2str(icell,'%02d-'),extractAfter(filenamesavepng,[stimPara.date,'\'])];
    else
        filename = rgcname('Reversing Grating WVSP', rgdata.clusters(icell,:),rgdata.savingpath);
        pngfilename = [num2str(icell,'%02d-'),filename];
    end
    p.title(filename);
    
    savepngFast(h,rgdata.savingpath, pngfilename);
    
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', icell,Ncells,toc);
    fprintf(msg);
end
close(h);
%----------------------------------------------------------------------

end
