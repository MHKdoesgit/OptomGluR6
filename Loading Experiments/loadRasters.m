

function newExperiment = loadRasters( experiment, expId, varargin )
%
%%% loadRasters %%%
%
%
% This function Loads and bin spikes from rasters files.
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
% Copied from same function from Fernando (2011-11-29),
% modified by Mohammad, 30.07.2014
% update to faster version by omiting data cell-structure on 14.06.2015


stimPara = experiment.stimPara;
clusters = experiment.clusters';
ftimes = experiment.ftimes;

numcells = size(clusters, 2);

if (isfield(stimPara,'nblinks'))
    nblinks = stimPara.nblinks;
else
    disp('There is no nblinks in stimPara file, binning using ftimes.');
    nblinks = ceil(mean(diff(ftimes))/(1/60)); %blinks for 60Hz Monitor in seconds
end

if isnan(nblinks) || nblinks==0, nblinks = 1; end   % to avoid error when there is no ftimes
numframes = nblinks * numel(ftimes);

spikesbin = zeros( (numframes/nblinks), numcells );
%spikesPerFrame = zeros( numframes, numcells );
spktimes = cell(1,size(clusters,2));

idx = 1;
for clus = clusters
    dataFile = sprintf('%d_SP_C%d%02d.txt', expId, clus(1), clus(2));
    try
        spiketimes = load([experiment.originalFolder '/rasters/' dataFile], '-ascii');
    catch ME
        disp(ME.message)
        warning(['file ',dataFile,'not found, brace yourself crash is coming!']);
    end
    spktimes{idx} = spiketimes;
    tmp = spikesBin(ftimes, spiketimes);
    spikesbin(:,idx) = tmp;
    idx = idx+1;
end

experiment.spiketimes = spktimes;
experiment.spikesbin = spikesbin;
experiment.nblinksBin = nblinks;

newExperiment = experiment;

end
