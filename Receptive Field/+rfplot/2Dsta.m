

function varargout = plot2Dsta(mainpath, varargin)
%
%%% plot2Dsta %%%
%
%
% This function make a surf plot of 2D STA from checker flicker stimulus.
% Becuase imagesc cannot show the polarity of the response of the cell but
% here by using surf function we can see the polarity and then location of
% the responce. This function also plot all 5 important plot of the 3
% different colors for comparison between the cells.
%
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   Plot (1): for all 20 frame of each color.
%   Plot (2): 5 important plot of all 3 colors.
%
%   Note : this function uses the RFData output of Calc_RF function and for
%   plotting and saving the plot it usese subplot_tight and savepng
%   internaly.
%
% written by Mohammad, 13.08.2015


totaltime =tic;
if (nargin < 1)
    mainpath = uigetdir();
end;

folderpath = [mainpath, '\Data Analysis\RF_Analysis\'];

foldername = dir(folderpath); foldername = {foldername(3:end).name};

if numel(foldername) > 1
    for ii=1:numel(foldername), disp(foldername{ii});end;
    prompt = 'Too many checkerflicker on the dance floor!! choose only one!! ==> ';
    foldername = foldername(input(prompt));
end;

folderpath = [folderpath,foldername{1},'\'];

datapath = [folderpath,'RF_data\'];
datafiles = dir([datapath,'*.mat']);    datafiles = {datafiles.name};


if ~exist([folderpath,'\RF_sta2dplot'],'dir')
    mkdir ([folderpath,'\RF_sta2dplot']);
    mkdir ([folderpath,'\RF_sta_all']);
end;

plotpos = {[1 12],[3 14],[5 16],[7 18],[9 20]};
plpos = [0,5,10];

for ii = 1:numel(datafiles)
    
    load([datapath,datafiles{ii}]);
    
    stimuliColors = fieldnames(RFData);
    if numel(stimuliColors) == 1
        colmp = {'redblue'};
    elseif numel(stimuliColors) == 3
        colmp = {'red', 'green', 'purple'};
    end;
    
    h1 = figure(1);
    for jj = 1:numel(stimuliColors)
        %%
        stimSize = size(RFData.(stimuliColors{jj}).STA);
        bestFrame = RFData.(stimuliColors{jj}).svd.max3ind;
        if bestFrame < 5, bestFrame = 5;end;
        if bestFrame > stimSize(3)-1, bestFrame = stimSize(3)-1;end;
        
        signalRange = ceil(max([max(max(RFData.(stimuliColors{jj}).STA(:,:,(bestFrame)))),...
            abs(min(min(RFData.(stimuliColors{jj}).STA(:,:,(bestFrame)))))])*10)/10;
        
        figure(1);
        for ww = 1:5
            subplot_tight(numel(stimuliColors),5,plpos(jj)+ww,0.04)
            surf(RFData.(stimuliColors{jj}).STA(:,:,(bestFrame-4)+ww),...
                'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
            %colormap(brighten(purple,-0.5))
            colormap(colmp{jj});
            caxis([-signalRange signalRange])
            shading flat;
            view(-40, 4);
            axis([0 stimSize(2) 0 stimSize(1) -signalRange signalRange]);
            axis square;
            title(['frame No.: ',num2str((bestFrame-4)+ww)],'FontSize',12);
            set(gca,'xTick',0:stimSize(2)/2:stimSize(2),'yTick',0:stimSize(1)/2:stimSize(1));
            set(gca,'zTick',-signalRange:signalRange/2:signalRange);
            freezeColors;
        end;
        
        h2 = figure(2);
        for ww = 1:5
            subplot_tight(4,10,plotpos{ww},0.02)
            surf(RFData.(stimuliColors{jj}).STA(:,:,(bestFrame-4)+ww),...
                'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
            %colormap(brighten(purple,-0.5))
            colormap(colmp{jj});
            caxis([-signalRange signalRange])
            shading flat;
            view(-40, 4);
            axis([0 stimSize(2) 0 stimSize(1) -signalRange signalRange]);
            axis square;
            title(['frame No.: ',num2str((bestFrame-4)+ww)],'FontSize',12);
            set(gca,'xTick',0:stimSize(2)/2:stimSize(2),'yTick',0:stimSize(1)/2:stimSize(1));
            set(gca,'zTick',-signalRange:signalRange/2:signalRange);
            freezeColors;
        end;
        
        
        for kk = 1: stimSize(:,3)
            subplot_tight(4,10,20+kk,0.02)
            
            surf(RFData.(stimuliColors{jj}).STA(:,:,kk),'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
            colormap(colmp{jj});
            shading flat;
            caxis([-signalRange signalRange])
            axis([0 stimSize(2) 0 stimSize(1) -signalRange signalRange]);
            axis square;
            title(['frame No.: ',num2str(kk)],'FontSize',12);
            set(gca,'xTick',0:stimSize(2)/2:stimSize(2),'yTick',0:stimSize(1)/2:stimSize(1));
            set(gca,'zTick',-signalRange:signalRange:signalRange);
            view(-40, 4);
            freezeColors;
        end;
        set(h2, 'position',[100 20 1800 850]);
        set(h2,'color',[1 1 1]);
        suptitle_mod(h2,['Spike Triggered Average of ',stimuliColors{jj},' Stimulus',datafiles{ii}(8:end-4)],10);
        picfr = getframe(h2);
        imwrite(picfr.cdata, [folderpath,'RF_sta2dplot\',['Spike Triggered Average of', datafiles{ii}(8:end-4)],'.tiff'], 'WriteMode', 'append');
        close(h2);
        clear picfr;
    end;
    
    set(h1, 'position',[100 50 1800 900]);
    set(h1,'color',[1 1 1]);
    suptitle_mod(h1,['Spike Triggered Average of', datafiles{ii}(8:end-4)],10);
    picfr = getframe(h1);
    savepng(picfr.cdata,[folderpath,'RF_sta_all\',['Spike Triggered Average of', datafiles{ii}(8:end-4)],'.png'],10,900);
    close (h1);
    clear picfr RFData bestFrame h1 h2 signalRange jj kk ww;
end;

disp(seconds2human (toc(totaltime)));

load gong.mat;
sound(y);