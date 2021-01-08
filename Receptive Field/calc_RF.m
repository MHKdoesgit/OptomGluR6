

function varargout = calc_RF(datapath, varargin)
%
%%% calc_RF %%%
%
%
% This function calculate spatio-temporal receptive field of ganglion cells
% based on thier reponse to checker flicker stimulus
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
% based on Calc_Receptive_field and MEA_Receptive_Field and multiple
% functions used by Daisuke to calculate receptive field,
% written by Mohammad, 14.08.2014
% update to complete new version based on setSVDoutput and fit_RF on
% 06.06.2016.
% update to new version with faster execution time by calculating all the
% STAs together and fit them in one run. Update on 24.01.2017.
% small update for colormaps and cut of the RF plots based on thier 3-5
% sigma of the fit, done on 24.01.2017.
% update to use calculateBlockSTAbw function to have much faster STA
% calculation speed compared to the mex function on 14.05.2018.
% moved function fitRFstas to a separate function on 09.11.2018.


totaltime =tic;
if (nargin < 1), datapath = uigetdir(); end    % for no input opens a uigetdir
if (nargin > 1), cell2analyze = varargin{1}; else,  cell2analyze = 0; end
% get all the STAs and other parameters
tic;
[RFdata,clusters,para,savingpath,nameset,lcol] = getRFparameters(datapath,1,cell2analyze);
toc;
% fit all the STAs in one run and then plot them later
%gaussFit = fitRFstas(STA,normSTA,numSpikes,para);
% set color map
cmpcol = {(flipud(cbrewer('div', 'RdBu',255))),(flipud(cbrewer('div', 'RdBu',255)))};
visplot = 'off';
%wb = waitbarTimed(0,['Plotting receptive field data for cell 1 out of ',int2str(size(clusters,1))]);

parfor ii = 1: size(clusters,1)
    warning('off','MATLAB:prnRenderer:opengl');
    % clearvars RFData STA normSTA numSpikes;
    figCounter = 4;
    % waitbarTimed(ii/size(clusters,1),wb,['Plotting receptive field data for cell ',int2str(ii),' out of ',int2str(size(clusters,1))]);
    
    for jj = 1:numel(nameset)
        
        % putting data into structure
        %   RFData.(nameset{jj}) = gaussFit{jj}(ii);
        %   RFData.(nameset{jj}).para = para;
        %%% Plotting STA
        plotSTA( RFdata(ii,jj), cmpcol, clusters(ii,:), para, savingpath, visplot, jj);
        %%% plotting spatial and temporal components
        figCounter = plotBestSC(RFdata(ii,jj),'full',para,cmpcol,lcol, figCounter, visplot,jj);
        figCounter = plotBestSC(RFdata(ii,jj),'zoom',para,cmpcol,lcol, figCounter, visplot,jj);
    end
    
    % saving final plot
    filename = generateRGCname('Receptive field',clusters(ii,:),savingpath,1);
    %     filename = {['Receptive field of cell ',num2str(clusters(ii,1)),', cluster ',num2str(clusters(ii,2))];[' for experiment on ',para.expdate]};
    joinRFplots(savingpath,filename,para,visplot);
    % saving the receptive field data
    % if strcmp(para.colormode,'monochromatic'), RFData = RFData.blackwhite; end
    % save([savingpath,'\rf_data\',filename{1},filename{2}],'-struct','RFData');
    % clearvars figCounter filename
    disp(['Analysis of channel: ',num2str(ii),' is finito']);
    
end
warning('on','MATLAB:prnRenderer:opengl');
%close(wb);
disp(' And BOOM!!! Goes all the Receptive Fields ');
disp(seconds2human (toc(totaltime)));
disp('Hallelujah hallelujah hallelujah hallelujah hallelujah');
sound(struct2array(load('handel.mat','y')))
varargout{1} = RFdata;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [gaussFit,clus,para,savingpath,ns,col,colmp] = getRFparameters(dp, saving,varargin)

if nargin> 2, num2analyze = varargin{1}; else, num2analyze = 0; end

expData = loadRawData(dp,'checkerflicker','RF_Analysis','rf_data','excludename','1600x1200'); % only load data, no folder making

if isempty((regexpi(dp,{'color','marmoset'}))), colormode = 'monochromatic';  else,  colormode = 'dichromatic';   end % checking color mode
if ~any(isfield(expData.stimPara,{'color','Color'})), colormode = 'monochromatic'; end
chsize = [num2str(expData.stimPara.stixelheight),'x',num2str(expData.stimPara.stixelwidth)];    % checking stixel size
if expData.stimPara.blackwhite, chType = '_binary'; else,    chType = '_gauss';  end    % checking stimuli type bi or gauss

[savingpath,para.seed,para.meanintensity,para.contrast,ns,col,colmp,~] = setWhiteNoiseParameters(expData.stimPara,[dp,'/Data Analysis/'],colormode);
savingpath = [savingpath,'_',chsize,chType,'_',num2str(expData.stimPara.nblinks),'blinks'];
para.colormode = colormode;

if saving
    %%%folder making
    if ~exist(savingpath,'dir'), mkdir(savingpath); end
    if ~exist([savingpath,'/rf_sta'],'dir')
        mkdir ([savingpath,'/rf_sta']);
        mkdir ([savingpath,'/rf_fit']);
        mkdir ([savingpath,'/rf_data']);
    end
end

% Stimuli Parameters
para.nblinks= expData.stimPara.nblinks;
para.blackwhite = expData.stimPara.blackwhite;
para.resolution = floor(1e3/60*para.nblinks);      % ==> up-sampling ratio, change this for different resolution.
timebehindspikeinSec = 680;   %  700 ms behind each spike
resbyFrame = round(floor(1e3/60*para.nblinks)/para.resolution);
para.timebins = floor(timebehindspikeinSec/(1e3/60*para.nblinks)) *resbyFrame;
if (para.timebins/resbyFrame) < 15, para.timebins = 15 * resbyFrame; end % minimum bins are 15
para.time = fliplr(linspace(0,(para.resolution*para.timebins), para.timebins));
para.screen = expData.screen;
para.stixelWidth = expData.stimPara.stixelwidth;
para.stixelHeight = expData.stimPara.stixelheight;
para.Nx = ceil(para.screen(1,1)/para.stixelWidth);
para.Ny = ceil(para.screen(1,2)/para.stixelHeight);
para.rfsigma = 1.5;
para.expdate =  expData.date;
switch expData.lightprojection
    case 'LightCrafter'
        para.monitorpixel = 8;
    case 'OLED'
        para.monitorpixel = 7.5;
end

% loading clusters
clus = expData.clusters;
if num2analyze(1) > 0  % this is to do partial analysis
    clus = clus(num2analyze(1):num2analyze(2),:);
end
% loading frametimes
ft = expData.ftimes;
if para.nblinks > 2     % in case nblinks were higher than 2
    ft = ft(1:para.nblinks/2:end);
end
% for up-sampling
if resbyFrame > 1
    ft = linspace(ft(1),ft(end),length(ft)*resbyFrame);
end
% loading spike times
% spkbin = histc(CelltoMatUE(expData.spiketimes(1:size(clus,1))),ft);
spk = expData.spiketimes;
spkbin = zeros(length(ft),size(clus,1));
for i=1:size(clus,1)
    if isempty(spk{i})
        spk{i} = 0;%normalizetoRange(rand(1,length(ft)),min(ft),max(ft));
    end
    spkbin(:,i) = histc(spk{i},ft);
end   % binnig all spikes
numspk = cellfun(@numel,expData.spiketimes(1:size(clus,1)));

% calculating STA for all the cells and all the colors
STA = calcRFsta(spkbin,para);

% fit all the STAs in one run and then plot them later
[RFdata, gaussFit] = fitRFstas(STA,numspk,para, ns);
% saving all the data here before going for plotting
saveRFdata(savingpath,clus,strcmp(para.colormode,'monochromatic'),RFdata);

end

%---------------------------------------------------------------------------------------------------%

function [out,badspk,fcout] = checkSpikesQuality(para, fc,ch,clus)  %#ok

if para.numSpike < 50       % when the spikes are fewer than 50 skip the calculation and use NaN
    warning([num2str(para.numSpike),' spikes?? What is it?!? Spikes for Ants??','No analysis for channel: ',num2str(ch),',cluster: ',num2str(clus)]);
    [out.STA,out.normFactor,out.tempComp, out.spatialComp,out.gaussfit,out.centercoord,out.correctedcenter, out.peakpos,...
        out.peakvalue, out.subrow,out.subcol,out.RFdiameter] = deal(nan(1,1));
    out.para = para;
    badspk = true;
    figure(fc); set(gcf,'visible','off','position',[400 100 610 935],'color',[1 1 1]); plot(NaN);
    figure(fc+1); set(gcf,'visible','off','position',[400 100 610 935],'color',[1 1 1]); plot(NaN);
    fcout = fc + 2; % little trick to avoid code crash in joining plot part.
else
    disp(['You have ', num2str(para.numSpike),' spikes, Be hold! the STA is coming!']);
    badspk = false;    out = NaN;    fcout = fc;
end

end

%---------------------------------------------------------------------------------------------------%

function STA3D = calcRFsta(spkbin,para, varargin)

% reshape function to make 3D STA
%stafun = @(x) (permute(reshape(x,para.timebins,para.Ny,para.Nx),[2 3 1]));
% Euclidean norm of STA
%normSTAfun = @(x) (sqrt(sum(bsxfun(@power,x,2),3)));   % for 2-D norm
STA3D  = cell(1,numel(para.seed));

% for ii = 1:numel(para.seed)
%     if para.nblinks == 2
%         rawSTA = CalcSTAHR(spkbin,para.timebins,para.nblinks,para.Nx*para.Ny,para.seed(ii),...
%             para.blackwhite,para.meanintensity(ii),para.contrast(ii),para.resolution);
%     elseif para.nblinks == 1
%         rawSTA = CalcSTAHR(spkbin,para.timebins,para.nblinks*2,para.Nx*para.Ny,para.seed(ii),...
%             para.blackwhite,para.meanintensity(ii),para.contrast(ii),para.resolution*2);
%     else            % to solve isues with higher nblinks (fix this later)
%         rawSTA = calcSTA_mex(spkbin,para.timebins,para.nblinks/2,para.Nx*para.Ny,para.seed(ii),...
%             para.blackwhite,para.meanintensity(ii),para.contrast(ii));
%     end
%     STA3D{ii} = cellfun(stafun,rawSTA,'UniformOutput',0);
%     normSTA{ii} = reshape(cell2mat(cellfun(normSTAfun,STA3D{ii},'UniformOutput',0)),para.Ny,para.Nx,size(spkbin,2));
% end

numcells = size(spkbin,2);
chunksize = factor(size(spkbin,1)); % the smaller the chunk the closer to real ouptut the STA is for cutting at the start of each chunk
chunksize = chunksize(ceil(end/2)); % one smaller than middle factor is selected for the chunk size
% for the old cpu function uncomment these lines to have better memory
% management!
% mem = memory;   % this is to make sure that the stimulus array will fit in the memory
% if (para.timebins*para.Nx*para.Ny*size(spkbin,1)/chunksize * 4) > mem.MaxPossibleArrayBytes
%     chunksize = factor(size(spkbin,1));
%     chunksize = chunksize(ceil(end/2)+1);
% end
if chunksize == size(spkbin,1)
    spkbinchunks = reshape(spkbin',numcells,chunksize,[]);
else
    spkbinchunks = reshape(spkbin',numcells,[],chunksize);
end
stafun = @(x)(reshape(x,para.Ny,para.Nx,para.timebins));

for ii = 1:numel(para.seed)
    
    rawSTA = calculateSTAbwGPU(spkbinchunks,para.timebins,para.Nx*para.Ny,para.seed(ii),para.contrast(ii));
    STA3D{ii} = cellfun(stafun,num2cell(double(rawSTA),[2 3]),'UniformOutput',0);
    % for the old cpu function uncomment these lines
    % rawSTA = calculateBlockSTAbw(spkbinchunks,para.timebins,para.Nx*para.Ny,para.seed(ii),para.contrast(ii));
    % STA3Dall = reshape(rawSTA',numcells, para.Ny, para.Nx,para.timebins);      % to have 4D sta use this combination
    % STA3D{ii} = cellfun(stafun,num2cell(rawSTA,1),'UniformOutput',0);
    
    %normSTA{ii} = reshape(cell2mat(cellfun(normSTAfun,STA3D{ii},'UniformOutput',0)),para.Ny,para.Nx,size(spkbin,2));
end

end

%--------------------------------------------------------------------------------------------------%

function plotSTA(rfdata, cmp, ch, para, path, vis, iter)
cmp = cmp{iter};
% disp('Plotting STA');
h1 = figure('position',[250 100 para.timebins*70 850],'color',[1 1 1],'visible',vis);
fs2D = round(max([abs(max(rfdata.STA(:))),abs(min(rfdata.STA(:)))])/0.1)*0.1;     % for z-axis
scSize = para.screen;

for kk = 1:para.timebins
    subplot_tight(4,para.timebins/4,kk,[0.01 0.02])
    imagesc(1:para.stixelWidth:scSize(1),1:para.stixelHeight:scSize(2),rfdata.STA(:,:,kk));
    plotRFcoordinates(rfdata.correctedcenter,'allcolors',rgb('gray'),'lineslinewidth',0.25,'circlelinewidth',0.5);
    axis equal tight off;    axis([0 scSize(1) 0 scSize(2)]); box off; colormap(cmp);  caxis([-fs2D fs2D]);
    %xticks(0:scSize(1)/2:scSize(1)); yticks(0:scSize(2)/2:scSize(2));
    %set(gca,'Xtick', 0:scSize(1)/2:scSize(1), 'YTick', 0:scSize(2)/2:scSize(2),'fontsize',7);
    %if (mod(kk,para.timebins/4)==1), axis on;  set(gca,'xcolor','w');    end
    %if (kk>para.timebins-(para.timebins/4)), axis on; set(gca,'ycolor','w'); end
    %if (kk==para.timebins-(para.timebins/4)+1), set(gca,'xcolor','k','ycolor','k'); end
    title(['STA frame No.',num2str(kk)],'fontSize',9);
    if kk==rfdata.peakpos(3), title(['Best STA frame No.',num2str(kk)],'fontSize',8,'color',rgb('orange'));end
end
filename = ['STA for cell ',num2str(ch(1)),', cluster ',num2str(ch(2)), ', for experiment on ', para.expdate];
suptitle_mod(h1,filename,4);
%saving the plot
if strcmpi(para.colormode,'monochromatic')
    savepngFast(h1,[path,'\rf_sta\'],filename);
else
    imwrite(print2array(h1, 2, '-opengl'), [path,'\rf_sta\',filename,'.tiff'], 'WriteMode', 'append','Compression','lzw');
end
close(h1);
end

%--------------------------------------------------------------------------------------------------%

function figCounter = plotBestSC(rfdata,plotType,para,cmp,col, figCounter, vis,iter)

cmp = cmp{iter};
col = col{iter};

%regvals = 2*pi*((1:1e3)-1)/(1e3-1);
%regfun = @(sig)(sig * sqrtm(rfdata.RFdiameter)* [cos(regvals);sin(regvals)] + repmat(rfdata.gaussfit.correctmu,size(regvals)));

h(figCounter) = figure(figCounter);
set((h(figCounter)),'visible',vis,'position',[400 50 610 935],'color',[1 1 1]);

switch lower(plotType)
    case {'full','all','non-zoom','big'}
        fullregion = floor((rfdata.RFdiameter*5)/(para.monitorpixel*para.stixelHeight)/2);
        if fullregion < 1, fullregion = 5; end  % to avoid curshes in case of noisy data
        x_img = rfdata.peakpos(2)-fullregion:rfdata.peakpos(2)+fullregion;
        y_img = rfdata.peakpos(1)-fullregion:rfdata.peakpos(1)+fullregion;
        x_img = x_img(x_img>=1 & x_img<=para.Nx);
        y_img = y_img(y_img>=1 & y_img<=para.Ny);
        SC = rfdata.STA(y_img,x_img,rfdata.peakpos(3));
        screenrow = 1:para.stixelHeight:para.screen(1);
        screencol = 1:para.stixelWidth:para.screen(2);
        x_img = screenrow(x_img);
        y_img = screencol(y_img);
        mi = maxIndex(SC);
        f2Drange = 0.01;
        if (SC(mi(1),mi(2)) < 0), SC = -1*SC; end
        if skewness(SC(:))<0, SC = -1*SC; end % spatial component should be always positive
        % sigline = regfun(5);
        
    case {'zoom','fit','small'}
        x_img = rfdata.subrow;
        y_img = rfdata.subcol;
        SC = rfdata.spatialComp;
        f2Drange = 0.1;
        % sigline = regfun(3);
end
%x_img = x_img(x_img > min(sigline(1,:)) & x_img < max(sigline(1,:)));
%y_img = y_img(y_img > min(sigline(2,:)) & y_img < max(sigline(2,:)));

subplot(15,10,[1 60]);
plot(para.time,rfdata.tempComp,'color',rgb(col),'LineWidth',2);
axis([0 600 -0.8 0.8]);    grid on; axis square;
xticks(0:200:800);      yticks(-1:0.5:1);       set(gca,'fontsize',8);
title ('Temporal Component','FontSize',10);
ylabel('Filter Strength','FontSize',8);

subplot(15,10,[72.5 88]);
sumX = sum(SC,1)/ max(abs(sum(SC,1)));
plot(linspace(x_img(1),x_img(end),length(sumX)),sumX,'o-','color',rgb(col),'lineWidth',1,'markerfacecolor',rgb(col),'markersize',4);
set(gca,'xaxisLocation','top');
try
    axis([x_img(1) x_img(end) floor(min(sumX)/0.2)*0.2 0.01+1.1*ceil(max(sumX)/0.2)*0.2]);
    xticks(ceil(x_img(1:5:end)/10)*10);      yticks(-2:0.5:2);
catch ME
    disp(ME.message)
    axis([0 1 -2 2]);
end
set(gca,'fontsize',8);
title (['Center diameter: ', num2str(round(rfdata.RFdiameter,1)),', Surround diameter: ',num2str(round(rfdata.surrounddiameter,1))] ,'FontSize',9);

subplot(15,10,[91 148]);
f2D = ceil(max([abs(max(SC(:))),abs(min(SC(:)))])/f2Drange)* f2Drange;
imagesc(x_img,y_img,SC);
clb = colorbar; colormap(gca,cmp);  caxis([-f2D f2D]); % freezeColors;
clb.Location = 'westoutside';
clbY = get(clb,'YTick');        set(clb,'YTick',min(clbY):mean(clbY)-min(clbY):max(clbY));
hold on;
plot(rfdata.correctedcenter(1,:),rfdata.correctedcenter(2,:),'color',rgb('black'),'LineWidth',1.5);
plotRFcoordinates(rfdata.correctedsurround,'allcolors',rgb('dimgray'),'lineslinewidth',0.5,'circlelinewidth',1);
%plot(sigline(1,:),sigline(2,:),'color',rgb('brown'),'LineWidth',1.5);
axis equal tight;     set(gca,'Xtick',[],'YTick',[]);

subplot(15,10,[99 150]);
sumY = sum(SC,2)/ max(abs(sum(SC,2)));
plot(sumY,linspace(y_img(1),y_img(end),length(sumY)),'o-','color',rgb(col),'lineWidth',1,'markerfacecolor',rgb(col),'markersize',4);
set(gca,'yaxisLocation','right');       grid on;
try
    axis([floor(min(sumY)/0.2)*0.2 0.01+1.1*ceil(max(sumY)/0.2)*0.2, y_img(1) y_img(end)]);
    xticks(-2:0.5:2);           yticks(ceil(y_img(1:5:end)/10)*10);
catch ME
    disp(ME.message);
    axis([-2 2 0 1]);
end
set(gca,'fontsize',8);
%
% arranging plots based on the main RF plot
% rfsp.Position = [rfsp.Position(1), rfsp.Position(2)-0.02, rfsp.Position(3),rfsp.Position(4)];
% tcsp.Position = [rfsp.Position(1), tcsp.Position(2), rfsp.Position(3),rfsp.Position(4)];
% rfup.Position = [rfsp.Position(1), rfsp.Position(2)+rfsp.Position(4)+0.01, rfsp.Position(3),rfup.Position(4)];
% rfdn.Position = [rfsp.Position(1)+rfsp.Position(3)+0.02, rfsp.Position(2), rfdn.Position(3),rfsp.Position(4)];
% clb.Location = 'westoutside';
% clb.Position = [clb.Position(1),rfsp.Position(2),clb.Position(3)-0.01,clb.Position(4)];
% clbY = get(clb,'YTick');        set(clb,'YTick',min(clbY):mean(clbY)-min(clbY):max(clbY));
% t.Position = [t.Position(1), t.Position(2)+0.05, t.Position(3)];
set(gcf,'position',[400 100 610 935]);
sortoftightfig(gcf,0.7);
figCounter = figCounter +1;

end

%--------------------------------------------------------------------------------------------------%

function joinRFplots(Spath,name,para,vis)

switch para.colormode
    case 'monochromatic'
        for kk = 4:5, sortoftightfig(figure(kk),0.5); end
        savefinalfigure(figure(5),Spath,name,[500 100 480 840],vis);
        savefinalfigure(figure(4),[Spath,'\rf_fit'],name,[500 100 480 840],vis);
    case 'dichromatic'
        for kk = 4:7, sortoftightfig(figure(kk),0.7); end
        h1 = combinefigs(figure(4),figure(6),'LR');
        savefinalfigure(h1,[Spath,'\rf_fit'],[name{1},name{2}],[100 10 1120 985],vis);
        h2 = combinefigs(figure(5),figure(7),'LR');
        savefinalfigure(h2,Spath,[name{1},name{2}],[100 10 1120 985],vis);
    case 'trichromatic'
        for kk = 4:9, sortoftightfig(figure(kk),0.7); end
        h10 = combinefigs(figure(4),figure(6),'LR');
        h1 = combinefigs(h10,figure(8),'LR');
        savefinalfigure(h1,[Spath,'\rf_fit'],[name{1},name{2}],[100 10 1120 985],vis);
        h20 = combinefigs(figure(5),figure(7),'LR');
        h2 = combinefigs(h20,figure(9),'LR');
        savefinalfigure(h2,Spath,[name{1},name{2}],[100 10 1120 985],vis);
end
%close all; % not needed because it closes the waitbar too.
end

%--------------------------------------------------------------------------------------------------%

function savefinalfigure(h, Spath,name, plotsize, vis)

set(h,'color',[1 1 1],'visible', vis);
suptitle_mod(h, name,3);
set(h,'position',plotsize);
if size(name,1)>1, name = [name{1},name{2}]; end
savepngFast(h,Spath,name);
%saveas2(h,[Spath,'\',name],300,'pdf');
%saveas2(h,[Spath,'\',name],300,'fig');
close(h)
end

%--------------------------------------------------------------------------------------------------%

function RFData = saveRFdata(dp,clus,monochromflag,expinput)

fns =fieldnames(expinput);
for ii = 1:size(clus,1)
    if monochromflag
        RFData = expinput.blackwhite(ii);
    else
        for jj = 1:numel(fns)
            RFData.(fns{jj}) = expinput.(fns{jj})(ii);
        end
    end
    filename = generateRGCname('Receptive field',clus(ii,:),dp);
    save([dp,'\rf_data\',filename],'-struct','RFData');
end
end
