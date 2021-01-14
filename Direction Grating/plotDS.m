

function [x,y,circAvg] = plotDS (ds, varargin)
%
%%% plotDS %%%
%
%
% This function plots the polar plot for showing the preferred direction of
% the response.
%
%================================Inputs====================================
%
%   colords : DS plot color.
%   colordsspot : DS spot colors.
%   colorsumvector : color of the vector sum.
%   colorsumvectorspot : color of the spot at the tip of the vector sum.
%   colorpatchedge : edge color of the patch plot.
%   linewidthds : line width of the DS plot.
%   linewidthsumvector : line width of the vector sum line.
%   sizedsspot : marker size of the DS plot spots.
%   sizewidthsumvector : marker size of the vector sum plot.
%   drawpatch :flag for drawing the patch.
%   drawdsspot : flag for drawing DS spots.
%   drawsumvectorspot : flag for drawing vector sum marker.
%   facealpha : transparency of the patch plot.
%   gridlines : line style of the grid.
%   gridlinewidth : line width of the grid.
%   gridcolor : color of the grid lines.
%   gridratio : ratio for grid line compared to the circles (default is 1).
%   fontsize : font size of the firing rate labels.
%
%================================Output====================================
%
%   DS plot : the main output is DS plot.
%   x,y,circAvg : the coordinate of the DS plot and the vector sum.
%
% written by Mohammad, 16.02.2018 based on the same function from Fernando
% that was used since 2015.


%defaults:____________________
def.colds = [197 53 24]/255; % stimulus duration in sec
def.coldsspot = [245 164 136]/255;    % pre-frame duration in sec
def.colsumvec = [91 137 199]/255; % green on, uv off color
def.colsumvecspot = [29 113 184]/255; % green off, uv on color
def.colpatchedge = [197 53 24]/255;
def.lwds = 1.5;   % pre-frame color
def.lwsumvec = 2;  % transparency
def.sizedsspot = 6;   % open new figure
def.sizevecsumspot = 8;  % this is to flip the oreder of the plot to match right-to-left data on integration plot
def.drawpatch = true;   % draw base line below each bar
def.drawdsspot = false;   % linewidth for all the stairs plots
def.drawvecsumspot = true;% width of the shaded region (in seconds after preframes)
def.facealpha = 0.5;
def.gridline = '-';
def.gridlw = 0.5;
def.gridcol = [1 1 1]*0.4;
def.gridratio = 1;
def.fontsize = 7.5;


% parsing inputs: _____________________
pnames = {'colords','colordsspot','colorsumvector','colorsumvectorspot','colorpatchedge','linewidthds','linewidthsumvector',...
    'sizedsspot','sizewidthsumvector', 'drawpatch','drawdsspot','drawsumvectorspot','facealpha',...
    'gridlines','gridlinewidth','gridcolor','gridratio','fontsize'};

[para.colds, para.coldsspot, para.colsumvec, para.colsumvecspot,para.colpatchedge, para.lwds, para.lwsumvec, ...
    para.sizedsspot, para.sizevecsumspot, para.drawpatch, para.drawdsspot, para.drawvecsumspot,...
    para.facealpha, para.gridline, para.gridlw, para.gridcol, para.gridratio, para.fontsize] = ...
    internal.stats.parseArgs(pnames, struct2cell(def)', varargin{:});

% calculate whatever is needed
angs = ds.anglesRep;
frate = ds.perAngleRep;% / sum(data.perAngleRep);

x = frate .* cos(angs);
y = frate .* sin(angs);

circAvg = circularAverage(ds.perAngle(:),ds.angles);
circAvg = circAvg*sum(ds.perAngle);

% Draw the custom polar plot grid
if all(frate<0)
    radiiresp = linspace(2*floor(min( ds.perAngle)/2),0,5 );
    radii = linspace(2*floor(min(ds.perAngle)/2),0, 5);
else
    radiiresp = linspace(0, 2*ceil(max( ds.perAngle)/2),5 );
    radii = linspace(0, 2*ceil(max(ds.perAngle)/2), 5);
end
%drawGrid( radii, numel(angs)-1 , para.gridline);
drawDSGrid( radii, numel(angs)-1,  para.gridline, para.gridcol, para.gridlw, para.gridratio);

% Draw axes
iter = length(radiiresp);
for r = radii(end:-2:2)
    text( r, 0, {sprintf(' %g', round(radiiresp(iter)));' Hz'},'verticalalignment', ...
        'middle','horizontalalignment', 'left','fontsize',para.fontsize);
    iter = iter-1;
end

hold on;
if para.drawpatch
    patchfun(x,y,para.colds,para.colpatchedge,para.facealpha, para.lwds);
else
    plot(x, y, '-','color', para.colds, 'linewidth', para.lwds,'markersize',para.sizedsspot);
end

if para.drawdsspot
    plot(x, y, 'o','markerfacecolor', para.coldsspot, 'markeredgecolor', 'none','linewidth', para.lwds,'markersize',para.sizedsspot);
end
plot([0 circAvg(1)], [0 circAvg(2)], 'color', para.colsumvec, 'linewidth', para.lwsumvec);
if para.drawvecsumspot
    plot(circAvg(1), circAvg(2), 'o','color', para.colsumvecspot,'markerfacecolor', para.colsumvec,...
        'markersize', para.sizevecsumspot,'linewidth', para.lwsumvec);
end
axis tight off;
hold off;


end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function p = patchfun(x,y,col,coledge,transp, lw)

if logical(transp)
    p = patch([x,fliplr(x)],[y,fliplr(y)],1,'facecolor',col,'edgecolor',coledge,'facealpha',transp,'linewidth',lw);
else
    p = patch([x,fliplr(x)],[y,fliplr(y)],1,'facecolor',col,'edgecolor',coledge,'linewidth',lw);
end

end

%---------------------------------------------------------------------------------------------------%

function drawDSGrid( radii, nangles, lstyl, lincolor, lw, gridratio, varargin )
radii = abs(radii); % to avoid cases with negative firing rate difference from background
for r = radii
    %if r <= 0,         continue;    end
    rectangle('Position', [-r, -r, 2*r, 2*r],'Curvature', [1 1], 'linestyle', lstyl,'edgecolor', lincolor,'linewidth',lw);
end

angles = linspace(0,360,nangles+1);
angles = angles(1:(end-1));

for ang = angles
    line([0 1]*gridratio*cosd(ang)*max(radii),[0 1]*gridratio*sind(ang)*max(radii),...
        'linestyle', lstyl, 'color', lincolor, 'linewidth',lw);
end
axis square;

end
