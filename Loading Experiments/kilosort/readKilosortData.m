

function ksrasters = readKilosortData(ksdir, varargin)
%
%%% readKilosortData %%%
%
%
% This function reads all the numpy arrays (*.NPY) files generated by
% kliosort and modifed by phy to .mat format. Additionally, this function
% separate the spike times and other data for each stimulus and channel.
% The output is ready to use as a normal raster file.
%
%================================Inputs====================================
%
%   ksdir : path to the kilosorted data folder. Must contain all the files.
%
%================================Output====================================
%
%   ksrasters : rasters per channels, per stimulus,
%
% written by Mohammad, 14.07.2019.
% Updated to match phy 2.0 output and exclude pca_features on 15.10.2019.
% re-wamped andupdate to new version for using with phy 2.0a on 19.02.2020.

tic;

% setting some options
p = inputParser();
p.addParameter('loadpcfeatures', false, @(x) islogical(x));
p.addParameter('matchtsvfiles', true, @(x) islogical(x));
p.addParameter('loadallnpys', false, @(x) islogical(x));
p.addParameter('alignalltemplates', 'good', @(x) ischar(lower(x)));
p.addParameter('saveextra', false, @(x) islogical(x));
p.parse(varargin{:});
ops = p.Results;

if ~contains(ksdir,'ks_sorted'), ksdir = [ksdir,'/ks_sorted/']; end

% make a folder for kilosort rasters
savingpath = [ksdir,'/ksrasters/'];
if ~exist(savingpath, 'dir'), mkdir(savingpath); end

% first loading the params file to get sorting parameter
ksraw.params_py = loadParamsPy(fullfile(ksdir, 'params.py'));

largenpyfiles = {'spike_clusters.npy','spike_templates.npy', 'spike_times.npy','amplitudes.npy',...
    'pc_features.npy','template_features.npy'};
smallnpyfiles = dir([ksdir,'/*.npy']);      smallnpyfiles = {smallnpyfiles.name};
smallnpyfiles = smallnpyfiles(not(ismember(smallnpyfiles, largenpyfiles)))';

% loading spike cluster. This file has the indices for all the clusters,
% the indices is used for separartion into clusters and stimuli.
spike_clusters = readNPY(fullfile(ksdir, 'spike_clusters.npy'));

% defining the function for separating spike vectors into clusters.
[ksraw.clusters_id,~,clusnumids] = unique(spike_clusters);
ksgroupfun = @(x)(accumarray(clusnumids,1:size(x,1),[],@(r){x(r,:)}));
% https://www.mathworks.com/matlabcentral/answers/345111-splitting-matrix-based-on-a-value-in-one-column

% loading and separating big files first.
if ops.loadallnpys
    ksraw.spike_clusters = ksgroupfun(spike_clusters);
end
%clearvars spike_clusters;   % clean up each for memory efficiency


spike_templates = readNPY(fullfile(ksdir, 'spike_templates.npy'));
ksraw.spike_templates = ksgroupfun(spike_templates);
%clearvars spike_templates;

% loading spike time in sampling rates
spike_times = readNPY(fullfile(ksdir, 'spike_times.npy'));
%ksraw.spike_times = ksgroupfun(spike_times);
%clearvars spike_times;
% loading spike amplitudes
amplitudes = readNPY(fullfile(ksdir, 'amplitudes.npy'));
ksraw.amplitudes = ksgroupfun(amplitudes);
%clearvars amplitudes;

% loading templates
%temps = readNPY(fullfile(ksdir, 'templates.npy'));     % loaded as part of small npys below

% Now loading smaller npy files
% npynamelist = {'channel_map.npy','channel_positions.npy','pc_feature_ind.npy','similar_templates.npy',...
%     'template_feature_ind.npy','templates.npy','templates_ind.npy','whitening_mat.npy','whitening_mat_inv.npy'};
for ii = 1:numel(smallnpyfiles)
    ksraw.(smallnpyfiles{ii}(1:end-4)) =  readNPY(fullfile(ksdir, smallnpyfiles{ii}));
end
ksraw.channel_map = ksraw.channel_map + 1;  % +1 for MATLAB indexing and to match with KSoutput


% loading comments here first to correct misaligned templates
switch lower(ops.alignalltemplates)
    case 'all'
        comments = 'all';
    case 'good'
        comments = readClusterGroupsCSV(fullfile(ksdir, 'cluster_group.tsv'));
        comments = comments(:,2);
    case 'misaligned'
        comments = readClusterGroupsCSV(fullfile(ksdir, 'cluster_comment.tsv'));
        comments = comments(:,2);
end

[spktimes, avgtemps] = fixMisalignedUnits(spike_times, clusnumids, spike_templates, ksraw.templates, comments);

ksraw.spike_times = ksgroupfun(spktimes);
ksraw.templates = avgtemps;
% cleaning all the variables
clearvars spike_templates spike_times amplitudes spktimes avgtemps comments;

if ops.loadpcfeatures
    % loading PCA features
    if exist(fullfile(ksdir, 'pc_features.npy'),'file')
        pc_features = readNPY(fullfile(ksdir, 'pc_features.npy'));
        ksraw.pc_features = ksgroupfun(pc_features);
        clearvars pc_features;
    end
    
    if exist(fullfile(ksdir, 'template_features.npy'),'file')
        template_features = readNPY(fullfile(ksdir, 'template_features.npy'));
        ksraw.template_features = ksgroupfun(template_features);
        clearvars template_features;
    end
end


% loading tsv files and read their content.
tsvnames = dir([ksdir,'/*.tsv']);
tsvnames = {tsvnames.name};

for ii =  1:numel(tsvnames)
    cvals = readClusterGroupsCSV(fullfile(ksdir, tsvnames{ii}));    % use the modified version from MHK
    ksraw.(tsvnames{ii}(1:end-4)) = cvals;
end
clearvars -except ksraw ksdir savingpath ops;
fprintf('loading npy, tsv files, elapsed time: %2.2f s...\n', toc);

if ~isfield(ksraw,'cluster_quality') && isfield(ksraw,'cluster_qual') % to match naming issues between phy 1 and phy 2
    ksraw.cluster_quality = ksraw.cluster_qual;
    ksraw = rmfield(ksraw,'cluster_qual');
end

if ops.matchtsvfiles  % this is to have same format as the cluster_info.tsv file
    
    if isfield( ksraw.cluster_info,'Amplitude') % for ks2
        ksraw.cluster_Amplitude = [ksraw.cluster_info.id, ksraw.cluster_info.Amplitude];
    end
    if isfield( ksraw.cluster_info,'ContamPct') % for ks2
        ksraw.cluster_ContamPct = [ksraw.cluster_info.id, ksraw.cluster_info.ContamPct];
    end
    if isfield( ksraw.cluster_info,'KSLabel')   % for ks2
        ksraw.cluster_KSLabel   = [num2cell(ksraw.cluster_info.id), ksraw.cluster_info.KSLabel];
    end
    if isfield( ksraw.cluster_info,'amp')       % for ks1
        ksraw.cluster_amp = [ksraw.cluster_info.id, ksraw.cluster_info.amp];
    end
    ksraw.cluster_comment   = [num2cell(ksraw.cluster_info.id), ksraw.cluster_info.comment];
    ksraw.cluster_group     = [num2cell(ksraw.cluster_info.id), ksraw.cluster_info.group];
    ksraw.cluster_quality   = [ksraw.cluster_info.id, ksraw.cluster_info.quality];
end

% this is to remove cells that are put to noise after checking the repsonse
% of the cells. all the qualities that have the noise label are put to 5.
noisyafteranalysis = (ksraw.cluster_quality(:,2) < 5 & contains(ksraw.cluster_group(:,2),'noise'));
ksraw.cluster_quality(noisyafteranalysis,2) = 5;
ksraw.cluster_info.quality(noisyafteranalysis) = 5;

% sort based on channel and cluster_ids
[chnums,chsorted] = sortrows([[ksraw.cluster_info.ch],[ksraw.cluster_info.id]],[1 2],'ascend');
clus = zeros(size(chnums,1),1);
for ii = 1:size(chnums,1)
    clusidx  = eq(chnums(:,1),ii);
    clus(chsorted(clusidx)) = 1:sum(clusidx);
end
% this is just to add a cluster number for each id and channel to match igor naming system.
ksraw.cluster_info.clus = clus;

% here is to get the starting and ending point of each experiment
if exist(fullfile(ksdir, 'bininfo.mat'),'file') % first try to get it from bininfo file
    bininfo = struct2array (load(fullfile(ksdir, 'bininfo.mat')));
    stimsamplessum = cumsum( bininfo.stimsamples);
    stimsamplerates = [[1;stimsamplessum(1:end-1)+1],stimsamplessum];
    ksraw.stim_start_end = stimsamplerates;
    ksraw.sampling_rate = bininfo.fs;
else % if it was not in the bin info then read it from mcd file, add msrd option later!
    expfolder = extractBefore(ksdir,'ks_sorted');
    mcdpath = [expfolder,'/'];
    if exist(mcdpath,'dir')
        [ksraw.stim_start_end,fs] = getMCDsamplingrates(mcdpath); % needs Neuroshare functions
    else
        [ksraw.stim_start_end,fs] = getMCDsamplingrates();    % if the path is wrong open the GUI
    end
    ksraw.sampling_rate = fs;
end

% these used in the end to separate the spikes into stimuli
numexperiments = size(ksraw.stim_start_end,1);
expduration = double(ksraw.stim_start_end);
samplerate = double(ksraw.params_py.sample_rate);
if ~eq(samplerate, ksraw.sampling_rate)
    error('DA-FAQ, the shit is wrong with the sampling rate, check the params.py file!');
end

% now re-organizing raw data into array-structure for easier indexing
clusnums = length(ksraw.cluster_info.id);
ks = cell(clusnums,1);
clusqual = ksraw.cluster_quality;
ksinfofn = fieldnames(ksraw.cluster_info); % field names of ksraw.cluster_info

for ii = 1:clusnums
    % this sorting is essential since accumarray change the order of
    % indices!
    [~,sortedspkidx] = sort(ksraw.spike_times{ii},'ascend');
    
    ks{ii}.cluster_id = single(ksraw.clusters_id(ii));
    ks{ii}.spike_times = ksraw.spike_times{ii}(sortedspkidx);% ./ ksraw.params_py.sample_rate;
    ks{ii}.amplitudes = ksraw.amplitudes{ii}(sortedspkidx);
    if ops.loadallnpys
        ks{ii}.spike_clusters = ksraw.spike_clusters{ii}(sortedspkidx); % all the orders are corrected!
        ks{ii}.spike_templates = ksraw.spike_templates{ii}(sortedspkidx);
    end
    spktemplates = ksraw.spike_templates{ii}(sortedspkidx);
    if isfield(ksraw,'pc_features')
        ks{ii}.pc_features = ksraw.pc_features{ii}(sortedspkidx);
    end
    if isfield(ksraw,'template_features')
        ks{ii}.template_features = ksraw.template_features{ii}(sortedspkidx);
    end
    
    [ti,tc] = unique(spktemplates);   % this is to get the channel number and relevant templates
    
    ks{ii}.pc_feature_ind = ksraw.pc_feature_ind(ti+1,:);   % +1 for MATLAB indexing
    %ks{ii}.similar_templates = squeeze(ksraw.templates(ii,:,:));     % ti for old code without template aliging
    ks{ii}.template_feature_ind = ksraw.template_feature_ind(ti+1,:) +1;
    ks{ii}.templates = squeeze (ksraw.templates(ii,:,:));
    %ks{ii}.templates_ind = ksraw.templates_ind(ti+1,:,:);
    
    [minamp,minloc] = min(min(ks{ii}.templates));
    
    %     tempchvar = squeeze(var(ks{ii}.templates,[],1));
    %     if size(tempchvar,2) < size(tempchvar,1), tempchvar = tempchvar'; end
    %     [maxamp,maxloc] = max(tempchvar,[],2);
    %     [~, ampidx] = sort(maxamp,'descend');     % sort templates to get the most common ones!
    ks{ii}.template_index_loc = ti;
    ks{ii}.template_ind_frq = tc;
    ks{ii}.template_amplitude = minamp;        % easy way to find best template
    ks{ii}.template_order = minloc;%maxloc(ampidx)-1;
    ks{ii}.template_ch = ksraw.channel_map(minloc); %mode(maxloc(ampidx))  % get the most common template (this is for phy1, for phy2 you can use cluster_info.ch)
    ks{ii}.channel_number =  ksraw.cluster_info.ch(ii);
    %if ~isequal(maxloc, minloc), disp(['variance-peak mismatch in template: ', num2str(ks{ii}.template_ch)]); end
    if ~isequal(ks{ii}.channel_number,ks{ii}.template_ch)
        %if size(tempchvar,1) > 1, champidx = maxloc(ampidx); else,  [~, champidx] = sort(tempchvar(),'descend'); end
        %if numel(champidx) > 5, champidx = ksraw.channel_map(champidx(1:5))'; else, champidx = ksraw.channel_map(champidx)'; end
        minamp = min(ks{ii}.templates);      [~, champidx] = sort(minamp,'ascend');  champidx = ksraw.channel_map(champidx(1:5))';
        disp(['channel number: ',num2str(ks{ii}.channel_number),', template found at: ',num2str(champidx)]);
        clearvars champidx;
    end
    
    ks{ii}.cluster_group = ksraw.cluster_group(ii,:);
    % this is for backward compatibility especially phy1 and old version of readClusterGroupsCSV
    clqu = find(clusqual(:,1)==ks{ii}.cluster_id);
    if isempty(clqu), clqu = [ks{ii}.cluster_id,5]; else, clqu = clusqual(clqu,:); end % fill up null/empty indices
    ks{ii}.cluster_quality = clqu;
    
    for jj = 1:numel(ksinfofn)
        ks{ii}.cluster_info.(ksinfofn{jj}) = ksraw.cluster_info.(ksinfofn{jj})(ii);
    end
    % this is cosmetic, to see fast what shit is inside each cluster and
    % also for printing in excel
    ks{ii}.channel_id = [ks{ii}.cluster_id, size(ks{ii}.spike_times,1), ks{ii}.cluster_group,ks{ii}.cluster_quality, ks{ii}.channel_number];
    
    clearvars tempchannel sortedspkidx ti tc clqu tempchvar maxamp maxloc sortamp ampidx spktemplates jj;
    
end

ks = cell2mat(ks);

ksinfo = rmfield(ksraw,{'spike_templates','spike_times','amplitudes','pc_feature_ind',...
    'template_feature_ind','templates','whitening_mat','whitening_mat_inv'});
if ops.saveextra
    % saving the raw data (it is heavy and slow, be patient my child!)
    save([savingpath,'\ksrawdata.mat'],'-v7.3','-struct','ksraw');
    save([savingpath,'\ksinfo.mat'],'-v7.3','-struct','ksinfo');
end
fprintf('creating ksraw file, elapsed time: %2.2f s...\n', toc);

% re-arranging the raw data into different channels
[~,chsorted] = sortrows([[ks.channel_number];[ks.cluster_id]]',[1 2],'ascend'); % sort based on channel and cluster_ids
ksch = ks(chsorted');

clearvars -except ks ksdir ksch numexperiments expduration samplerate savingpath ksinfo ops;
fprintf('re-arranging data per stimulus & per channels, elapsed time: %2.2f s...\n', toc);
% separating each channel to its stimuli components

ftosplt = {'amplitudes','spike_times'};
if ops.loadpcfeatures,  ftosplt = [ftosplt,'spike_clusters','spike_templates'];    end
if isfield(ksch(1),'pc_features'),    ftosplt = [ftosplt,'pc_features'];    end
if isfield(ksch(1),'template_features'),    ftosplt = [ftosplt,'template_features'];    end

for kk = 1:numel(ftosplt)
    ksperexp = cell(size(ksch,1),numexperiments);
    thisfield = ftosplt{kk};
    
    for ii=1:size(ksch,1)
        
        for jj =1 : numexperiments
            spktime= (ksch(ii).spike_times);
            spkidx = spktime > expduration(jj,1) & spktime <= expduration(jj,2);
            
            if strcmpi(thisfield, 'spike_times')    % spike-times is changed to ms and subtracted from start of stimulus
                ksperexp{ii,jj} = (double(ksch(ii).(thisfield)(spkidx)) ./ samplerate) - (expduration(jj,1) ./ samplerate);
            else
                ksperexp{ii,jj} = ksch(ii).(thisfield)(spkidx);
            end
        end
        ksch(ii).(thisfield) = ksperexp(ii,:);
    end
    ksrasters.(thisfield) = ksperexp;
end

for ii=1:size(ksch,1)
    ksrasters.template_info(ii,1) = rmfield(ksch(ii),ftosplt);
    ksrasters.sort_info(ii,1) = ksch(ii).cluster_info;
    ksrasters.clusters(ii,:) = [ksch(ii).channel_number, ksch(ii).cluster_info.clus, ...
        ksch(ii).cluster_quality(:,2),ksch(ii).cluster_quality(:,1)];
end
fn =fieldnames(ksinfo);
ksrasters.sort_params = rmfield(ksinfo,fn(contains(fn,'cluster_')));
ksrasters.template_info = rmfield(ksrasters.template_info,'cluster_info');
ksrasters.sort_info = rmfield(ksrasters.sort_info,{'sh'});

% saving the final data into two formats
if ops.saveextra
    save([savingpath,'\ksrastersperchannel.mat'],'-v7.3','ksch');
end
save([savingpath,'\ksrasters.mat'],'-v7.3','-struct','ksrasters');

fprintf('finito...good luck with the analysis...total time ===> %2.2f mins...\n', toc/60);


end

