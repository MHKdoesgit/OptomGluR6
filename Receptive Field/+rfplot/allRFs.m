

function allRFs(allrfcoords, targetnum, scsize, varargin)
%
%%% plotallRFs %%%
%
%
% This function plot all the receptive fields of a recording together and
% highlight one of the selected ganglion cells. It also has the option to
% plot the spatial component of the receptive field in background of the
% gaussian fit.
%
%
% ===============================Inputs====================================
%
%   allrfcoords : receptive field coordinates of all the ganglion cells.
%   targetnum: target number of the ganglion cell.
%   scsize : screen size (should match the unit of the coordinates).
%   rfimage : receptive field image (optional).
%
%================================Output====================================
%
%   no output : this is a plotting function and no output is produced.
%
% written by Mohammad, 19.10.2018.
% update for new marmoset dataset on 03.03.2020.


% setting some options
p = inputParser();
p.addParameter('rfimage',[]);
p.addParameter('color', [0.25*[1 1 1],0.5]);
p.addParameter('targetcolor', [1 0 0]);
p.addParameter('lw', 0.5,@(x) isnumeric(x));
p.addParameter('targetlw', 1.5, @(x) isnumeric(x));
p.addParameter('outline', false, @(x) islogical(x));
p.addParameter('colormap', flipud(cbrewer('div','RdBu',255)));
p.addParameter('outlinecolor', 0.8*[1 1 1]);
p.parse(varargin{:});
pltops = p.Results;


if not(isempty(pltops.rfimage))
    imagesc(1:scsize(1)/size(pltops.rfimage,2):scsize(1),1:scsize(2)/size(pltops.rfimage,1):scsize(2),pltops.rfimage);
    f2D = max(pltops.rfimage,[],'all');
    axis equal tight off;
    caxis([-f2D f2D]);
    colormap(gca, pltops.colormap);
end

hold on;
if pltops.outline
    rectangle('pos',[0 0 scsize],'facecolor','none','edgecolor',pltops.outlinecolor,'linewidth',pltops.lw);
end

if iscell(allrfcoords)
    for ii = 1:size(allrfcoords,1)
        line(allrfcoords{ii}(1,:),allrfcoords{ii}(2,:),'color',pltops.color,'linewidth',pltops.lw);
    end
    line(allrfcoords{targetnum}(1,:),allrfcoords{targetnum}(2,:),'color',pltops.targetcolor,'linewidth',pltops.targetlw);
else
    line(squeeze(allrfcoords(1,:,:)),squeeze(allrfcoords(2,:,:)),'color',pltops.color,'linewidth',pltops.lw);
    line(allrfcoords(1,:,targetnum),allrfcoords(2,:,targetnum),'color',pltops.targetcolor,'linewidth',pltops.targetlw);
end

if nargin > 3
    axis tight off;
    set(gca,'ydir','normal','tickdir','out');
    %set(gca,'ydir','reverse','tickdir','out');
else
    axis equal on;
    axis([0 scsize(1) 0 scsize(2)]);
    xticks(0:scsize(1)/2:scsize(1));
    yticks(0:scsize(2)/2:scsize(2));
    set(gca,'ydir','normal','tickdir','out');
    %set(gca,'ydir','reverse','tickdir','out');
end

end
