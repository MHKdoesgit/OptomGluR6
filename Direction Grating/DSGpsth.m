

function dgout = DSGpsth(ft, spktimes, para, stimID, stimphasedelay,varargin)

% stimulus parameter
nangles = para.Nangles(stimID);
cycles = para.cycles(stimID);
period = para.period(stimID);
nrep = (para.duration(stimID) ./ para.period(stimID))-1; % first trial is not counted
binlength = 10/1e3; % 10 ms psth bin
binnum = floor((period/ para.fps)/binlength);
binvec = linspace(0, (period/ para.fps), binnum+1);

ftimes = ft(1:floor(length(ft)/(nangles * cycles))* (nangles * cycles));
ftimes = reshape( ftimes, [], nangles, cycles);

cycras = deal(cell(cycles,1));
rasters = cell(nrep, nangles);
trialpsth = zeros(nrep, nangles, binnum);
for np = 1: nrep
    
    t0 = squeeze ( ftimes (np,:,:)) - stimphasedelay{stimID}; % Ignores the first period of each angle
    % the first pulse already starts after the end of the first period!
    tf = (squeeze ( ftimes (np, : ,:)) + period/ para.fps) - stimphasedelay{stimID};
    
    for ang = 1:nangles
        for cycl = 1:cycles
            cycras{cycl} = spktimes( and(spktimes > t0(ang, cycl) , spktimes < tf(ang, cycl)))-t0(ang, cycl);
        end
        rasters{np,ang} = CelltoMatUE(cycras);
        if isempty(rasters{np,ang})
            p = zeros(1,binnum);
        else
            if size( rasters{np,ang},2) < 2
                p = (histc(rasters{np,ang},binvec)/size(rasters{np,ang},2))*(1/binlength);  % to fix situation with 1 spike
            else
                p = mean(histc(rasters{np,ang},binvec),2)*(1/binlength);
                % 1-sum((p(:,1:2:end)- p(:,2:2:end)).^2,'all')./sum((p(:,1:2:end)-repmat(mean(p(:,1:2:end),2),[1 size(p(:,1:2:end),2)])).^2,'all');
            end
        end
        trialpsth(np,ang,:) = p(1:binnum);   % must transpose contspk for histc to work properly.
    end
end

% ratesOdd = squeeze(psth(1:2:end,1,:));
% ratesEven = squeeze(psth(2:2:end,1,:)); ratesOdd = ratesOdd(1:size(ratesEven,1),:);
% qualityRsq = squeeze(1-sum((ratesOdd-ratesEven).^2,3)./sum((ratesOdd - repmat((mean(ratesOdd,3)),1,1,size(ratesOdd,3))).^2,3));
%  qualityRsq = 1-sum((ratesOdd-ratesEven).^2,2)./sum((ratesOdd-repmat(mean(ratesOdd,2),[1 size(ratesOdd,2)])).^2,2);

resppsth = squeeze(mean(trialpsth,1));

Y = fft(resppsth);
P2 = abs(Y/binnum);
P1 = P2(:,1:floor(binnum/2)+1);
P1 = 2*P1;
f = (1/binlength)*(0:(binnum/2))/binnum;
f2 = para.fps/(binnum); f1=f2/2;
[~,f1ind] = min(abs(f-f1));
[~,f2ind] = min(abs(f-f2));

allF2 = P1(:,f2ind);
allF1 = P1(:,f1ind);

angles = rad2deg (rem(pi+(0:(nangles-1)) * 2*pi/nangles,2*pi));
[~,si] = sort(rad2deg(angles));
allF2 = allF2(si);
allF1 = allF1(si);
%angles = angles(si);

% allRatio = allF2./allF1;

dgout.trialrasters = rasters;
dgout.psth = resppsth;
dgout.trialpsth = trialpsth;
dgout.F1 = allF1;
dgout.F2 = allF2;
% dgout.angles = angles;

end
