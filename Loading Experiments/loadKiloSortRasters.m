
function newExperiment = loadKiloSortRasters( experiment, expId, ksrasters, varargin )
%
%%% loadKiloSortRasters %%%
%
%
% This function Loads and bin spikes from Kilosort rasters files.
%
% ===============================Inputs====================================
%
%    experiment : path to experiment folder
%    expId : Experiment number which used to load it's spike times
%
%================================Output====================================
%
%   newExperiment : structure containing spiketimes, binned spiketimes and
%                   nblinks of the experiment file.
%
% written by Mohammad, 18.06.2014 based on loadRasters function.
% updated for KS2 and phy2 on 19.02.2020.


if nargin < 3 || isempty(ksrasters)
    ksrasters = struct2array(load([experiment.originalFolder,'/ksrasters/ksrasters.mat']));
    warning('Yo, ksraster input is missing, you are loading all the raster for each experiment, this is slow as shit!');
end

if nargin > 3, biningflag = varargin{1}; else, biningflag = false; end % options to bin spikes with ftimes.

stimPara = experiment.stimPara;
clusters = experiment.clusters';
ftimes = experiment.ftimes;

numcells = size(clusters, 2);
if (isfield(stimPara,'nblinks'))
    nblinks = stimPara.nblinks;
else
    disp('There is no nblinks in stimPara file, binning using ftimes.');
    nblinks = round(mean(diff(ftimes))/(1/60)); %blinks for 60Hz Monitor in seconds
end

if isnan(nblinks), nblinks = 1; end   % to avoid error when there is no ftimes

ksclusters = ksrasters.clusters;
selectedclusters = ismember( ksclusters(:,3),experiment.clusters(:,3)); % use cluster ids for comparison
% to avoid index mixing between matlab and python phy.

% load all the rasters in once
spktimes = ksrasters.spike_times(selectedclusters,expId);

if biningflag
    numframes = nblinks * numel(ftimes);
    spikesbin = zeros( (numframes/nblinks), numcells );
    % binning like the old time, unecessary but kept for legacy codes
    for ii = 1: size(spikesbin,2)
        %tmp = spikesBinner_V2(ftimes, spktimes{ii});
        %tmpPF = spikesBinPerFrame(ftimes, spktimes{ii}, nblinks);
        spikesbin(:,ii) = spikesBin(ftimes, spktimes{ii});
    end
end

experiment.spiketimes = spktimes';
if biningflag
    experiment.spikesbin = spikesbin;
    % experiment.spikesPerFrame = spikesPerFrame;       % not used
    experiment.nblinksBin = nblinks;
end
experiment.amplitudes = ksrasters.amplitudes(selectedclusters,expId)';  % amplitudes from KS2
experiment.sortinginfo = ksrasters.sort_info(selectedclusters,:);    % cluster_info from phy2
experiment.sortinginfo = rmfield(experiment.sortinginfo,{'amp'});
newExperiment = experiment;

end
