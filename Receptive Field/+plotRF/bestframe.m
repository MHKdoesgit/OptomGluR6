

function bestframe(STA, peakpos, correctedcenter, correctedsurround, mor_I, para, pltops)
%
%%% plotRF.bestframe %%%
%
%
% This function plot the best frame of the sta along with receptive field
% fit of the center and surrond. It also make two plots for the sums of the
% best frame along the x and y axis.
%
%
% ===============================Inputs====================================
%
%   STA : whole sta (20-40 frames). If STA is a structure all the values
%         are taken from the subfield of the structure.
%   peakpos : frame number corresponding to the peak of sta.
%   correctedcenter : coordinate of corrected RF center.
%   correctedsurround : coordinate of corrected RF surround.
%   mor_I : moran's I measured during RF calculation, used for defining cell
%           quality. Based on Jian's subunit paper threshold is at 0.25.
%           changed to 0.2 to better match the data!
%   para : para file from the receptive field fits.
%   options : for defining properites of the output plot, check below.
%
%================================Output====================================
%
%   no output : this is a plotting function and no output is produced.
%
% written by Mohammad, 03.03.2020.


if nargin == 1  && isstruct(STA)    
    bestfr = STA.STA(:,:,STA.peakpos(3));
    centcoord = STA.correctedcenter;
    surrcoord = STA.correctedsurround;
    p = STA.para;
    mori = STA.moransI;
    
elseif nargin == 2  && isstruct(STA) && iscell(peakpos)
    bestfr = STA.STA(:,:,STA.peakpos(3));
    centcoord = STA.correctedcenter;
    surrcoord = STA.correctedsurround;
    p = STA.para;
    mori = STA.moransI;
    pltops = peakpos;    
    
elseif nargin == 2  && isstruct(STA) && numel(peakpos)==1    
    iter = peakpos;
    bestfr = STA.STA{iter}(:,:,STA.peakpos(iter,3)); % get the best frame of the STA
    centcoord = squeeze(STA.correctedcenter(:,:,iter));
    surrcoord = squeeze(STA.correctedsurround(:,:,iter));
    p = STA.para;
    mori = STA.moransI(iter,:);
    
elseif nargin == 3  && isstruct(STA) && numel(peakpos)==1 && iscell(correctedcenter)
    iter = peakpos;
    bestfr = STA.STA{iter}(:,:,STA.peakpos(iter,3)); % get the best frame of the STA
    centcoord = squeeze(STA.correctedcenter(:,:,iter));
    surrcoord = squeeze(STA.correctedsurround(:,:,iter));
    p = STA.para;
    mori = STA.moransI(iter,:);
    pltops = correctedcenter;    
else    
    bestfr = STA(:,:,peakpos(3));
    centcoord = correctedcenter;
    surrcoord = correctedsurround;
    p = para;
    mori = mor_I;
end


% setting some options
pplt = inputParser();
pplt.addParameter('color', [0.58 0 0.83]);
pplt.addParameter('centercolor', [0.2 0.2 0.2]);
pplt.addParameter('linecolor', [0 0.75 1]);
pplt.addParameter('lw', 0.5,@(x) isnumeric(x));
pplt.addParameter('centerlw', 1.5, @(x) isnumeric(x));
pplt.addParameter('colorbar', true, @(x) islogical(x));
pplt.addParameter('colorbarsouth', false, @(x) islogical(x));
pplt.addParameter('colorbarpos', [0.2 0.03]);
pplt.addParameter('colormap', flipud(cbrewer('div','RdBu',255)));
pplt.addParameter('goodrfthresh', 0.2,@(x) isnumeric(x));
pplt.addParameter('plotratio', 0.15,@(x) isnumeric(x));
pplt.addParameter('plotgapscale', 8,@(x) isnumeric(x));
pplt.parse(pltops{:});
pltops = pplt.Results;



% threshold 0.25 from Jians paper
if length(mori)>1
    mris = abs(diff(sort(mori,'descend')));
    if (mris(1)-mris(2)>0.25), goodrf = true; else, goodrf = false; end
else
    if (mori(1)>0.25), goodrf = true; else, goodrf = false; end
end

f2D = ceil(max([abs(max(bestfr(:))),abs(min(bestfr(:)))])/0.01)* 0.01;
%mxrf = bestfr(STA.peakpos(iter,1),STA.peakpos(iter,2));
%if mxrf < 0, offrf = true; else, offrf = false; end 
% sumX = sum(abs(bestfr),1)./size(bestfr,1);
% sumY = sum(abs(bestfr),2)./size(bestfr,2);
sumX = sum(bestfr,1)./size(bestfr,1);
sumY = sum(bestfr,2)./size(bestfr,2);
maxsum = max(abs([sumX(:);sumY(:)]),[],'all');
sumX = double(sumX ./ maxsum);
sumY = double(sumY ./ maxsum);
%sumX = sumX - median(sumX);
%sumY = sumY - median(sumY);

if isfield(p,'stixelHeight')
    xsx = 1:p.stixelHeight:p.screen(1);
    xsy = 1:p.stixelWidth:p.screen(2);
else
    xsx = 1:p.stixelheight:p.screen(1);
    xsy = 1:p.stixelwidth:p.screen(2);
end

if goodrf
    % fitting sumX and sumY with diff of gaussians!
    % p= [mu STDcenter scaleSTDsurr ampcenter ampsurr]
    predFun = @(p,x) (p(4)*(exp(-(x-p(1)).^2 /(2*pi*p(2).^2)))/(2*pi*p(2).^2) - ...
        p(5)*(exp(-(x-p(1)).^2/(2*pi*(p(3)*p(2)).^2)))/(2*pi*(p(3)*p(2)).^2)); % diff of Gaussians
    optimopts= optimset('Display','off','TolX', 1e3*eps,'TolFun', 1e3*eps,'MaxFunEvals',2e4,'MaxIter',2e4);
    %if offrf, sumX = -sumX; sumY = -sumY; end
    
    yy = sumX - median(sumX);
    [ampx, ampxloc] = max(sumX);
    gx = [xsx(ampxloc) 50 2 ampx 0.01];
    paramsSumX = lsqcurvefit(predFun, gx, xsx,yy, [0 0 1 0 0], [1000 100 50 2 2], optimopts);
    %plot(xx,sumX,xx,predFun(paramsSumX,xx))
    
    yy = (sumY - median(sumY))';
    [ampy, ampyloc] = max(abs(sumY));
    gy = [xsy(ampyloc) 50 2 ampy 0.01];
    paramsSumY = lsqcurvefit(predFun, gy, xsy,yy, [0 0 1 0 0], [1000 100 50 2 2], optimopts);
    %plot(xsy,sumY,xsy,predFun(paramsSumY,xsy))
    %if offrf, sumX = -sumX; sumY = -sumY; end
end

imagesc(xsx,xsy,bestfr);
hold on;
plotRF.coordinates(surrcoord,'allcolors',pltops.color,'lineslinewidth',pltops.lw,'circlelinewidth',pltops.lw+0.5);
plot(centcoord(1,:),centcoord(2,:),'color',pltops.centercolor,'LineWidth',pltops.centerlw);

xreg = ((pltops.plotratio/pltops.plotgapscale)*p.screen(1)) + [p.screen(1), p.screen(1) + (p.screen(1)*pltops.plotratio)];
yreg = ((pltops.plotratio/pltops.plotgapscale)*p.screen(2)) + [p.screen(2), p.screen(2) + (p.screen(2)*pltops.plotratio)];

if goodrf
    xplt = normalizetoRange([sumX;predFun(paramsSumX,xsx)],yreg(1),yreg(2));
    yplt = normalizetoRange([sumY';predFun(paramsSumY,xsy)],xreg(1),xreg(2));
else
    xplt = normalizetoRange(sumX,yreg(1),yreg(2));
    yplt = normalizetoRange(sumY',xreg(1),xreg(2));
end

% sumX plot
plot(xsx,xplt(1,:),'color',pltops.linecolor,'LineWidth',pltops.lw);
if goodrf
   plot(xsx,xplt(2,:),'color',pltops.color,'LineWidth',pltops.lw); 
end
% sumY plot
plot(yplt(1,:),xsy,'color',pltops.linecolor,'LineWidth',pltops.lw);
if goodrf
   plot(yplt(2,:),xsy,'color',pltops.color,'LineWidth',pltops.lw); 
end

% % sumX plot
% sc = ceil((p.screen(2)/5)/50)*50; % scaling of line plots
% if not(goodrf), sc = sc/2; end
% plot(xsx, p.screen(2)+50+sc*sumX,'color',pltops.linecolor,'LineWidth',pltops.lw);
% if goodrf
%     plot(xsx, p.screen(2)+50+sc*predFun(paramsSumX,xsx),'color',pltops.color,'LineWidth',pltops.lw);
% end
% % sumY plot
% plot(p.screen(1)+50+sc*sumY,xsy,'color',pltops.linecolor,'LineWidth',pltops.lw);
% if goodrf
%     plot(p.screen(1)+50+sc*predFun(paramsSumY,xsy),xsy,'color',pltops.color,'LineWidth',pltops.lw);
% end

axis equal tight on;
xmaxlim = 10 + ceil(xreg(2)/5)*5; %+ceil(max(p.screen(1)+50+sc*sumY)/15)*15;
ymaxlim = 10 + ceil(yreg(2)/5)*5; %+ceil(max(p.screen(2)+50+sc*sumX)/20)*20;
axratio = max(xsx)/xmaxlim;
axis([0 xmaxlim 0 ymaxlim]);
rfax = gca;
rfax.TickDir = 'out';
rfax.XAxisLocation = 'top';
rfax.YAxisLocation = 'right';
box off;
xticks(0:p.screen(1)/2:p.screen(1));
yticks(0:p.screen(2)/2:p.screen(2));

colormap(gca,pltops.colormap);   caxis([-f2D f2D]); % freezeColors;

if pltops.colorbar
    clb = colorbar;
    if pltops.colorbarsouth
        clb.Location = 'southoutside';
    else
        clb.Location = 'southoutside';
        clb.Position = [rfax.Position(1), rfax.Position(2)-pltops.colorbarpos(1), rfax.Position(3)*axratio, pltops.colorbarpos(2)];
    end
    clbY = get(clb,'YTick');        set(clb,'YTick',min(clbY):mean(clbY)-min(clbY):max(clbY));
    clb.TickDirection = 'out';
    % clb.Box = 'off';
    %clb.Position(4) = clb.Position(4)*0.85;
end
end