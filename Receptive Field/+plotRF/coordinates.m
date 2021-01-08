

function coordinates(rfcor, varargin)
%
%%% plotRF.coordinates %%%
%
%
% This function plots the receptive field coordinates with two lines to the
% edges of the figure.
%
%================================Inputs====================================
%
%   rfcor : receptive field coordinates
%   circlecolor : ellipse color.
%   linecolor : plot line colors.
%   allcolors : one color to rule them of all.
%   circlelinewidth : ellipse line width.
%   lineslinewidth : linewidth of the side lines.
%   linestyle : line style of the side lines
%   fontsize : text font size.
%   normalyaxis : to flip or not flip the y-axis if imagesc
%   textcolor : text color.
%
%================================Output====================================
%
%   RF ellipse plot : the ellipse would be plotted on top of the rf center.
%
% written by Mohammad, 12.02.2018.

%defaults:____________________
def.colcircle = 'k';
def.collines = 'k';
def.allcolor = [];
def.lwcircle = 2;
def.lwlines = 0.5;
def.lscircle = '-';
def.lslines = ':';
def.fontsize = 7.5;
def.normalyax = 1;
def.textcol = 'k';
def.axisoff = true;

% parsing inputs: _____________________
pnames = {'circlecolor','linecolor','allcolors','circlelinewidth','lineslinewidth','circlelinestyle','linestyle',...
    'fontsize','normalyaxis','textcolor','axisoff'};

[para.colcircle, para.collines,para.allcolor, para.lwcircle, para.lwlines, para.lscircle, para.lslines,para.fontsize,...
    para.normalyax, para.textcol,para.axisoff] = internal.stats.parseArgs(pnames, struct2cell(def)', varargin{:});

if ~isempty(para.allcolor)
    para.colcircle = para.allcolor; 
    para.collines = para.allcolor;
    para.textcol = para.allcolor;
end

if not(ishold), hold on; end

xl = get(gca,'xlim');
yl = get(gca,'ylim');

plot(rfcor(1,:),rfcor(2,:),'color',para.colcircle,'linewidth',para.lwcircle,'linestyle',para.lscircle);
plot([ xl(2)  max(rfcor(1,:))],[median(rfcor(2,:)) median(rfcor(2,:))],'color',para.collines,...
    'linewidth',para.lwlines,'linestyle',para.lslines);
plot([median(rfcor(1,:)) median(rfcor(1,:))],[yl(1)  min(rfcor(2,:))],'color',para.collines,...
    'linewidth',para.lwlines,'linestyle',para.lslines);
plot([median(rfcor(1,:)) median(rfcor(1,:))],[yl(2)  max(rfcor(2,:))],'color',para.collines,...
    'linewidth',para.lwlines,'linestyle',para.lslines);
plot([xl(1)  min(rfcor(1,:))],[median(rfcor(2,:)) median(rfcor(2,:))],'color',para.collines,...
    'linewidth',para.lwlines,'linestyle',para.lslines);

txtpos = min([diff(xl)/20, diff(yl)/20]);
text(xl(1)+txtpos, double(median(rfcor(2,:))+ txtpos/2), num2str(round(median(rfcor(2,:)))),...
    'fontsize',para.fontsize,'vert','baseline','color',para.textcol);
text( double(median(rfcor(1,:)))+txtpos/2,yl(2)-txtpos, num2str(round(median(rfcor(1,:)))),...
    'fontsize',para.fontsize,'color',para.textcol);

if para.normalyax
    set(gca,'ydir','normal');
end

if para.axisoff
    axis off;
end

end