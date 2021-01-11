

function staframes(sta, nframes, dims, rfpara, varargin)
%
%%% plotRF.staframes %%%
%
%
% This function plot defined numbers of STA frames. It can rearange the
% selected frame in user-defined dimension (rows and columns).
%
%
% ===============================Inputs====================================
%
%   sta : whole sta (20-40 frames).
%   nframes : number of frames to be plotted from the end.
%   dims : [nrows, ncols] for the final plot.
%   rfpara : para file from the receptive field fits.
%   options : for defining properites of the output plot, check below.
%
%================================Output====================================
%
%   no output : this is a plotting function and no output is produced.
%
% written by Mohammad, 03.03.2020.

% setting some options
p = inputParser();
p.addParameter('gap', [1 1], @(x) isnumeric(x));
p.addParameter('stapeak', nan(1,3));
p.addParameter('outline', true, @(x) islogical(x));
p.addParameter('colormap', flipud(cbrewer('div','RdBu',255)));
p.addParameter('outlinecolor', 'k');
p.addParameter('peakframeoutlinecolor', 'r');
p.addParameter('peakframefontcolor', 'r');
p.addParameter('showframenumber', true, @(x) islogical(x));
p.addParameter('fontcolor', 'k');
p.parse(varargin{:});
pltops = p.Results;

if numel(pltops.gap) == 1, pltops.gap = [pltops.gap, pltops.gap]; end
% sta = g.STA{ii};

mx = max(abs(sta),[],'all');
nrows = dims(1);
if numel(dims)==1
    ncols = ceil(nframes / nrows);
else
    ncols = dims(2);
end

xl = linspace(0,rfpara.screen(1),rfpara.Nx);
yl = linspace(0,rfpara.screen(2),rfpara.Ny);
% get the indices correct, if the nrows and ncols don't match it will crash
% here
if isfield(rfpara,'timebins')
    rfloc = reshape(rfpara.timebins-nframes+1:rfpara.timebins,nrows,ncols);
else
    rfloc = reshape(rfpara.nonlinBinN-nframes+1:rfpara.nonlinBinN,nrows,ncols);
end

idx = 0;
yidx = 0;

for kk = 1:nframes
    xm = idx*max(xl); if idx~=0, xm = xm+idx*(20*pltops.gap(1)); end
    ym = yidx*max(yl); if yidx~=0, ym = ym+yidx*(20*pltops.gap(2)); end
    if mod(kk,ncols)==0, yidx = yidx+1;end
    staimg = uint8(255*((sta(:,:,rfloc(kk))/mx)+1)/2);
    imagesc( xm + xl,ym + yl,staimg);
    hold on;
    if pltops.outline
        if ~isnan(pltops.stapeak(3)) && pltops.stapeak(3)==rfloc(kk)
            reccol = pltops.peakframeoutlinecolor;
        else
            reccol = pltops.outlinecolor;
        end
        rectangle('pos',[xm ym max(xl) max(yl)],'cur',[0 0],'facecolor','none','edgecolor',reccol,'linewidth',0.025);
    end
    idx = idx+1;
    if idx == (nframes/nrows); idx = 0; end
    if pltops.showframenumber
        if ~isnan(pltops.stapeak(3)) && pltops.stapeak(3)==rfloc(kk)
            txtcol = pltops.peakframefontcolor;
        else
            txtcol = pltops.fontcolor;
        end
        text(xm + xl(end-2),ym + yl(2),num2str(rfloc(kk)),'color',txtcol,'VerticalAlignment','top',...
            'HorizontalAlignment','right','FontSize',7,'FontAngle','italic');
    end
end
axis([0 xm+xl(end) 0 ym+yl(end)]);
axis equal tight off;

caxis([0 255]);
%caxis([-mx mx]);
colormap(gca,pltops.colormap);

end