

function drawGrid( radii, nangles, varargin )


if nargin < 2,        nangles = 8;     end;
if nargin < 1,        radii = [0.5 1];  end;
if nargin > 2,  lstyl = varargin{1}; else   lstyl = '--';    end;
if nargin > 3,  handle = varargin{2}; else handle = gca;    end;

lincolor = [1 1 1] * 0.4;
if nargin > 3,
    axis(handle);
    cla;
end;
axis square;

for r = radii
    if r <= 0,         continue;    end;
    rectangle('Position', [-r, -r, 2*r, 2*r],'Curvature', [1 1], 'linestyle', lstyl,'edgecolor', lincolor,'linewidth',1);
end;

angles = linspace(0,360,nangles+1);
angles = angles(1:(end-1));

for ang = angles
    line([0 1]*1.2*cosd(ang)*max(radii),[0 1]*1.2*sind(ang)*max(radii),'linestyle', lstyl, 'color', lincolor, 'linewidth',1);
end;

end