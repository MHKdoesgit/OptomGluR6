

function varargout = Analyze_DirectionGratingSequence(datapath, varargin)
%
%%% Analyze_DirectionGratingSequence %%%
%
%
% This function analyze the direction grating stimulus. The idea is to
% calculate the spiking frequency during each period of each angle. That is
% the amount of time that is take for the first bar to cross the screen.
% Based on this one can see whether the ganglion cells spikes to one
% direction more than other.
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%
%================================Output====================================
%
%   DSGCData : cell structure containing the measured properties of direction
%            grating stimulus.
%   plot : 1 to 4 ploar plot depicting the direction preference of each
%          ganglion cell in response to different color.
%
%   Note : this function use plotDS, drawgrid inside.
%
% written by Mohammad, 22.10.2015 based on analyzedirectionGrating written
% by Fernando for MEA course. plotDS and draw also borrowed from Fernando.
% update to new version with cleaner naming and structure but no change to
% the main calculation, on 30.01.2017.
% major update by re-writing the main code and adding options for p-value
% calculation and better plotting on 06.02.2018.
% added the feature to subtract the baseline from the stimulus in
% calculation of the dsi and osi, added on 24.01.2019.
% minor update for the opto project on 13.01.2021.

totaltime =tic;
if nargin < 1,    datapath = uigetdir();      end
if nargin > 1, bkgsub = varargin{1}; else, bkgsub = false; end 

% first calculate everything in one go and then plot
[dsgcall, clusters, plotoptions] = DS_analysis_parameters(datapath, bkgsub);
plotDSdata(dsgcall,clusters,plotoptions{:});

sound(struct2array(load('gong.mat','y')))
disp(seconds2human (toc(totaltime)));
varargout{1} = dsgcall;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [dsgcdata, clus, plotoptions] = DS_analysis_parameters(dp, bkgsub, varargin)

if bkgsub % in case of background subtraction
    dsfn = 'DirectionGratingSequence_Analysis_(pfr_sub_250ms)';
else
    dsfn = 'DirectionGratingSequence_Analysis';
end

[DSGCexp, savingpath] = loadRawData(dp,'directiongratingsequence',dsfn);

para = DSGCexp.stimPara;
if isfield(para,'nangles'), para.Nangles = para.nangles; end
if isfield(para,'preframeStimuli'), para.preframestimuli = para.preframeStimuli;  end
if isfield(para,'squareWave'), para.squarewave = para.squareWave;  end
para = rmfield(para,{'nangles', 'preframeStimuli', 'squareWave'});
if ~isfield(para,'numstimuli'), para.numstimuli = numel(para.Nangles); end
if length(para.Nangles) < para.numstimuli
    para.Nangles = repmat(para.Nangles,1,para.numstimuli);  % for older experiments with one number for angle numbers
end
para.fps = DSGCexp.info.screen.refreshrate;
if isfield(para,'lmargin')
para = rmfield(para,{'lmargin','rmargin','bmargin','tmargin'});
end
para.thresholds.ds = [0.15 0.05];
para.thresholds.os = [0.15 0.05];
para.thresholds.quality = 0.5;
if para.numstimuli > 2
    para.thresholds.stimperiod = [100 12];
else
    para.thresholds.stimperiod = para.period;
end
if ~isfield(para,'regeneration'), para.regeneration = para.preframeAngles; end
para.date = DSGCexp.date;

% first period doesn't have pulse
stimcutft = cumsum([0,(para.duration./para.period).*(para.Nangles.* para.cycles) - (para.Nangles .* para.cycles)]);  

clus = DSGCexp(1).clusters;

% this is fix a bug for cases when the preframe is not divisible to period.
stimphasedelay = delayofphaseshiftbug(para);
% calculate all the ds responses in one go!
[dsgcdata, dsgcpsth] = deal(cell(size(clus,1), para.numstimuli));
dstic = tic;
for ii = 1:para.numstimuli
    eachstimft = DSGCexp.ftimes(1+stimcutft(ii):stimcutft(ii+1));
    for jj = 1:size(clus,1)
        dsgcdata{jj,ii} = analyze_DS_stimulus(eachstimft, DSGCexp.spiketimes{jj}, para, ii, stimphasedelay, bkgsub);
        dsgcpsth{jj,ii} = DSGpsth(eachstimft, DSGCexp.spiketimes{jj}, para, ii, stimphasedelay);
    end
end
disp(seconds2human(toc(dstic)));
dsgcdata = cell2mat(dsgcdata);
dsgcpsth = cell2mat(dsgcpsth);
dsquality = DSOS_candidates(dsgcdata, para);

dsosdata = rearrnage_dsdata(dsgcdata);
dsosdata = struct2struct(dsosdata, dsquality);
disp(clus(dsquality.dscandicates,:));
% writing down the candidates for DS and OS
dscand = clus(dsquality.dscandicates,1:4);
fileID = fopen([savingpath,'ds_os_candidates.txt'],'W');
fprintf(fileID,'\n\nDS candidates, thresholds: dsi: %0.2f, dsi_pval: %0.2f, response quality: %0.2f\n\n',...
    [dsquality.thresholds.ds,dsquality.thresholds.quality]);
for jj = 1:sum(dsquality.dscandicates)
   fprintf(fileID,'\t%1.0f\t%1.0f\t%1.0f\t%1.0f\n',dscand(jj,:));
end
fprintf(fileID,'\n\n\n\nOS candidates, thresholds: osi: %0.2f, osi_pval: %0.2f, response quality: %0.2f\n\n',...
    [dsquality.thresholds.os,dsquality.thresholds.quality]);
oscand = clus(dsquality.oscandicates,:);
for jj = 1:sum(dsquality.oscandicates)
   fprintf(fileID,'\t%1.0f\t%1.0f\t%1.0f\t%1.0f\n',oscand(jj,:));
end
fclose(fileID);

% this is for adding commets from kilosort to the plots
if isfield(DSGCexp, 'sortinginfo'),   sortinfo = DSGCexp.sortinginfo;  else,    sortinfo = [];  end
if para.numstimuli == 6, fh = 1050; elseif para.numstimuli==2; fh = 715; else, fh = 850; end
% options for plotting, for more detail check plotdsdata function
plotoptions = {'colords',rgb('azure'),'colordslines',rgb('deepskyblue'),'position',[10 10 1850 fh],...
    'savefigure',true,'savedata',false,'savealltogether',false,'figuresavepath',savingpath,...
    'datasavepath',savingpath,'figurname','Direction Grating Sequence',...
    'visible','off','sortinfo',sortinfo};

[~,allsavename] = fileparts(savingpath(1:end-1));
allsavename = strrep([allsavename,'_for_experiment_on_',para.date],'Analysis','analysis');
ds.dsgcdata = dsgcdata;
ds.dsgcpsth = dsgcpsth;
ds.dsquality = dsquality;
ds.dsosdata = dsosdata;
save([savingpath,allsavename],'-struct','ds','-v7.3');

end

%--------------------------------------------------------------------------------------------------%

function dgout = analyze_DS_stimulus(ft, spktimes, para, stimID, stimphasedelay, bkgsub,varargin)

% stimulus parameter
nangles = para.Nangles(stimID);
cycles = para.cycles(stimID);
period = para.period(stimID);
if para.regeneration(stimID)==0
    bkgdur = 0;
else
    bkgdur = para.regeneration(stimID)/4/para.fps; % the last quarter of the gray frames are considered
    if bkgdur <= 0.1, bkgdur = 0.25; end  % cannot be shorter that  250 ms
end

ftimes = ft(1:floor(length(ft)/(nangles * cycles))* (nangles * cycles));
ftimes = reshape( ftimes, [], nangles, cycles);
t0 = squeeze ( ftimes (1,:,:)) - stimphasedelay{stimID}; % Ignores the first period of each angle
% the first pulse already starts after the end of the first period!
tf = (squeeze ( ftimes (end, : ,:)) + period/ para.fps) - stimphasedelay{stimID};

% pre-allocate variables
[nSpikes, fRate, pfrSpikes, pfrRate, firsttrialSpikes, firsttrialRate] = deal(zeros ( nangles, cycles));
[cycras,pfrcycras, firsttrialcycras, r] = deal(cell(cycles,1));
[rasters, pfrrasters, firsttrialrasters, rall] = deal(cell(nangles,1));


for ang = 1:nangles
    for cycl = 1:cycles
        cycras{cycl} = spktimes( and(spktimes > t0(ang, cycl) , spktimes < tf(ang, cycl)))-t0(ang, cycl);
        
        % to get the preframe, we should go back 1 period behind the first pulse and get preframe for duration of one period
        % this is because the pulse comes after the first period,so we should go two period behind.
        stimstart = t0(ang, cycl) - (period/para.fps);
        pfrdur = [stimstart - bkgdur, stimstart ];
        firsttrialdur = [stimstart , t0(ang, cycl) ];
%         if para.period(stimID) > para.regeneration(stimID) % in case the gray period before is shorter than one period
%             pfrdur = [stimstart - (para.regeneration(stimID)/para.fps), stimstart ];
%         else
%             pfrdur = [stimstart - para.period(stimID)/para.fps, stimstart ];    % if period is shorter than regeneration
%         end
        pfrSpikes( ang, cycl) = sum( and(spktimes > pfrdur(1) , spktimes < pfrdur(2)));
        pfrRate( ang , cycl ) = pfrSpikes( ang, cycl) ./ (diff(pfrdur));
        pfrcycras{cycl} = spktimes( and(spktimes > pfrdur(1) , spktimes < pfrdur(2)))-pfrdur(1);
        
        nSpikes( ang, cycl) = sum( and(spktimes > t0(ang, cycl) , spktimes < tf(ang, cycl)));
        fRate( ang , cycl ) = nSpikes( ang, cycl) ./ (tf ( ang, cycl )- t0( ang , cycl ));
        %  fRate( ang , cycl ) = fRate( ang , cycl ) - pfrRate( ang , cycl );  %subtraction from background
        
        % measuring same stuff for the first trial
        firsttrialSpikes( ang, cycl) = sum( and(spktimes > firsttrialdur(1) , spktimes < firsttrialdur(2)));
        firsttrialRate( ang , cycl ) = pfrSpikes( ang, cycl) ./ (diff(firsttrialdur));
        firsttrialcycras{cycl} = spktimes( and(spktimes > firsttrialdur(1) , spktimes < firsttrialdur(2)))-firsttrialdur(1);
        
        r{cycl} = spktimes( and(spktimes > pfrdur(1) , spktimes < tf(ang, cycl)))-pfrdur(1);
        
    end
    rasters{ang} = CelltoMatUE(cycras);
    pfrrasters{ang} = CelltoMatUE(pfrcycras);
    firsttrialrasters{ang} = CelltoMatUE(firsttrialcycras);
    rall{ang} = CelltoMatUE(r);
end

if bkgsub % for background subtraction
    fRate = fRate - mean(pfrRate(:)); %subtraction from average of all backgrounds
end
perAngle = squeeze( mean( fRate, 2 ) ); % Averages over the cycles
angles = rem(pi+(0:(nangles-1)) * 2*pi/nangles,2*pi); % ds grating seq starts from 180 degree
%(0:(nangles-1)) * 2*pi/nangles;
angles = angles'; % turns it into a column-vector
[~,angcorrectidx] = sort(rad2deg(angles));
angles = angles(angcorrectidx,:); % rearrange angles
perAngle = perAngle(angcorrectidx,:);

perAngleRep = [perAngle; perAngle(1, :)]; % Repeats the last angle to close the cycle
anglesRep = [angles; angles(1)];

circAvg = circularAverage(perAngle(:),angles);

dgout.totalSpikes = nSpikes(angcorrectidx,:);
dgout.fRate = fRate(angcorrectidx,:);
dgout.rasters = rasters(angcorrectidx,:);
% for preframe
dgout.pfrtotalSpikes = pfrSpikes(angcorrectidx,:);
dgout.pfrRate = pfrRate(angcorrectidx,:);
dgout.pfrrasters = pfrrasters(angcorrectidx,:);
% for the first trial
dgout.firsttrialSpikes = firsttrialSpikes(angcorrectidx,:);
dgout.firsttrialRate = firsttrialRate(angcorrectidx,:);
dgout.firsttrialrasters = firsttrialrasters(angcorrectidx,:);

dgout.rall = rall(angcorrectidx,:);

dgout.perAngle = perAngle;
dgout.angles = angles;
dgout.perAngleRep = perAngleRep;
dgout.anglesRep = anglesRep;
dgout.circAvg = circAvg*sum(perAngle);
%[dgout.dsi, dgout.osi] = calcDSIandOSI(dgout);
[dgout.dsi, dgout.dsi_angle, dgout.dsi_pval, dgout.osi, dgout.osi_angle, dgout.osi_pval,...
    dgout.dsi_dist, dgout.osi_dist] = getDirIndices(fRate(angcorrectidx,:)', angles', 1e3);  % function provided from Dimos
if dgout.dsi > 1    % to avoid weird dsi and osi when the firing rate is too low
    [dgout.dsi, dgout.dsi_angle] = deal(NaN);
end
if dgout.osi > 1
    [dgout.osi, dgout.osi_angle] = deal(NaN);
end
% separating parameters of each stimulus
dgout.para.stimulusID = stimID;
dgout.para.nangles = nangles;
dgout.para.cycles = cycles;
dgout.para.period = period;
dgout.para.duration = para.duration(stimID);
dgout.para.squarewave = para.squarewave(stimID);
dgout.para.gratingwidthwhite = para.gratingwidthwhite(stimID);
dgout.para.gratingwidthblack = para.gratingwidthblack(stimID);
dgout.para.regeneration = para.regeneration(stimID);
dgout.para.preframestimuli = para.preframestimuli;
dgout.para.contrast = para.contrasts(stimID);
dgout.para.numstimuli = para.numstimuli;
dgout.para.date = para.date;
dgout.para.fps = para.fps;

end

%--------------------------------------------------------------------------------------------------%

% function [dsi, osi] = calcDSIandOSI(dgin)
%
% dsi = abs(sum(exp(1i.* dgin.angles) .* dgin.perAngle/ sum( dgin.perAngle)));
% osi = abs(sum(exp((1i*2).* dgin.angles) .* dgin.perAngle/ sum( dgin.perAngle)));
%
% end

%--------------------------------------------------------------------------------------------------%

function dsos = DSOS_candidates(dsgcdata, para)

allCounts = reshape([dsgcdata.totalSpikes],size(dsgcdata(1).totalSpikes,1),...
    size(dsgcdata(1).totalSpikes,2),size(dsgcdata,1),para.numstimuli);
allCounts = permute(allCounts,[3 2 1 4]);

dsos.properstims    = para.period <= para.thresholds.stimperiod(1) & para.period >= para.thresholds.stimperiod(2); 
qualCounts     = allCounts(:, :, :, dsos.properstims);

dsos.responsequality = var(mean(qualCounts(:, :, :), 2), [], 3)./mean(var(qualCounts(:,:,:),[],3),2);

allDSIs = reshape([dsgcdata.dsi], size(dsgcdata,1), para.numstimuli);
pDSIs = reshape([dsgcdata.dsi_pval], size(dsgcdata,1), para.numstimuli);
allOSIs = reshape([dsgcdata.osi], size(dsgcdata,1), para.numstimuli);
pOSIs = reshape([dsgcdata.osi_pval], size(dsgcdata,1), para.numstimuli);

dsos.dscandicates  = sum(allDSIs(:,dsos.properstims) > para.thresholds.ds(1) ...
    & pDSIs(:,dsos.properstims) < para.thresholds.ds(2), 2)>1;
dsos.dscandicates(dsos.responsequality < para.thresholds.quality) = 0;

dsos.oscandicates  = sum(allOSIs(:,dsos.properstims) > para.thresholds.os(1) ...
    & pOSIs(:,dsos.properstims) < para.thresholds.os(2), 2)>1;
dsos.oscandicates(dsos.responsequality < para.thresholds.quality) = 0;
dsos.thresholds = para.thresholds;

end

%--------------------------------------------------------------------------------------------------%

function out = rearrnage_dsdata(dsgcdata)

[ncells, nstim] = size(dsgcdata);
[nang, ncyc] = size(dsgcdata(1).totalSpikes);

cf = @(x)(permute(reshape(x,nang,ncyc,ncells, nstim),[3 2 1 4]));
cf1d =  @(x)(permute(reshape(x,[],ncells, nstim),[2 1 3]));
cf1val = @(x) (reshape(x,ncells,nstim));

out.totalSpikes = cf([dsgcdata.totalSpikes]);
out.fRate = cf([dsgcdata.fRate]);
out.rasters = cf1d([dsgcdata.rasters]);
out.pfrtotalSpikes = cf([dsgcdata.pfrtotalSpikes]);
out.pfrRate = cf([dsgcdata.pfrRate]);
out.pfrrasters = cf1d([dsgcdata.pfrrasters]);
out.firsttrialSpikes = cf([dsgcdata.firsttrialSpikes]);
out.firsttrialRate = cf([dsgcdata.firsttrialRate]);
out.firsttrialrasters = cf1d([dsgcdata.firsttrialrasters]);
out.perAngle = cf1d([dsgcdata.perAngle]);
out.angles = cf1d([dsgcdata.angles]);
out.perAngleRep = cf1d([dsgcdata.perAngleRep]);
out.anglesRep = cf1d([dsgcdata.anglesRep]);
out.circAvg = cf1d([dsgcdata.circAvg]);
out.dsi = cf1val([dsgcdata.dsi]);
out.dsi_angle = cf1val([dsgcdata.dsi_angle]);
out.dsi_pval = cf1val([dsgcdata.dsi_pval]);
out.osi = cf1val([dsgcdata.osi]);
out.osi_angle = cf1val([dsgcdata.osi_angle]);
out.osi_pval = cf1val([dsgcdata.osi_pval]);
out.dsi_dist = cf1d([dsgcdata.dsi_dist]);
out.osi_dist = cf1d([dsgcdata.osi_dist]);
out.para = dsgcdata(1).para;
end
