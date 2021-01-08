

function [ thisexp ] = loadExperimentData( datapath , varargin)
%
%%% loadExperimentData %%%
%
%
% This function loads everthing about on expriment folder, this include all
% parameters of stimulus, goodcells and clusters, frametimings, spiketimings
% and binned spikes.
%
% ===============================Inputs====================================
%
%    datapath : path to experiment folder
%
%================================Output====================================
%
%   exp : multi-dimentional structure containing parameters of stimulus,
%         goodcells and clusters, frametimings, spiketimings
%         and binned spikes.
%
% based on loadExperiment function from Fernando,
% modified with new functions by Mohammad, 30.07.2014
% update to make indiviual exp files per experiment on 14.06.2015
% update the saving function for structures on 09.03.2016
% added new options for loading clusters with lower scores on 02.01.2017.
% added new options to load KiloSort data on 18.06.2019
% added new options for KS2 and phy2 on 19.02.2020.
% update to fully new version with input from the gui, more optimized data 
% structure, and new stimuli reading function etc on 08.01.2021.

if nargin < 1
    expinfo = loadExpInformation([]);
else
    if nargin > 1
        expinfo = loadExpInformation(datapath, varargin{1});
    else
        expinfo = loadExpInformation(datapath);
    end
end

% this part is to load correct rasters files
ksrastflag = false; igorrasflag = false;
if exist([datapath,'/ksrasters'],'dir'),    ksrastflag = true;  end
if exist([datapath,'/rasters'],'dir'),    igorrasflag = true;  end

% this is to select either kilosort or igor as you source of rasters
if ksrastflag && igorrasflag
    rastersource = questdlg('Choose your rasters','Kilosort or Igor','KiloSort','Igor','KiloSort');
    switch rastersource
        case 'KiloSort'
            ksrastflag = true; igorrasflag = false;
        case 'Igor'
            ksrastflag = false; igorrasflag = true;
        otherwise
            error('There aint no rasters here?!? Then whats the point of this shit!');
    end
end

if ksrastflag
    ksdata = load([datapath,'/ksrasters/ksrasters.mat']);
end

experiment.originalfolder = datapath;
[~, expname] = fileparts(experiment.originalfolder);
experiment.expname = expname;
experiment.date = expinfo.expdate;

% making folder for saving data
if not(exist(fullfile(datapath,'Data Analysis','Raw Data'),'dir'))
    mkdir(fullfile(datapath,'Data Analysis','Raw Data'));
end

% check frametimes path
if not(exist([datapath,'/frametimes/'],'dir'))
    ftpath = [datapath,'/ftimes/'];
else
    ftpath = [datapath,'/frametimes/'];
end
if not(exist(ftpath,'dir'))
    error('Yo, there aint no frametimes folder in this path, we cannt do shit without frametime!');
end

paramFile = [datapath '/parameters.txt'];
if not(exist(paramFile, 'file')) && expinfo.options.parameterstxtfile
    [files,header] = parameterstxtfromftimes(ftpath(1:end-1),expinfo.array.numelectrodes,expinfo.samplingrate);
    disp('loading experiment parameters');
    experiment.parameters.header = header;
    experiment.parameters.files = files;
    %experiment.parameters.numFiles = num;
end

disp('loading good cells');
experiment.clusters = loadClusters( datapath ,expinfo.spikesorting.bestscore, expinfo.spikesorting.worstscore);

if strcmpi(expinfo.spikesorting.datatype,'mcd') && strcmpi(expinfo.spikesorting.sorter,'igorpro')
    disp('get MEA coordinates');
    try
        meacoord = getMEAcoordinatesfromMCD(datapath);
        meaflag= true;
    catch ME
        disp(ME.message)
        warning('There aint no MEA pictures and MEA coordinates files here!');
        meaflag = false;
    end
else
    meaflag = false;
end


scres = expinfo.screen.resolution;

for ii = 1:numel(expinfo.stimulusnames)
    
    thisexp = experiment;
    
    disp('loading stimuli parameters');
    if strcmpi (expinfo.stimparatype,'auto')
        thisexp.stimPara = loadStimulusParameters(expinfo.stimulusnames{ii}, scres(1), scres(2));
    else
        try
            stimpath = [datapath,filesep,'stimuli',filesep,expinfo.stimulusnames{ii},'.txt'];
            thisexp.stimPara = readStimulusParameters(stimpath, scres(1), scres(2));
        catch ME
            disp(ME.message)
            warning('Da Fuuck! the stimulus is not registered to the readStimulusParameters function!');
            thisexp.stimPara = loadStimulusParameters(expinfo.stimulusnames{ii}, scres(1), scres(2));
        end
    end
    
    disp('loading frametimes'); % can be done cleaner but necessary to have 1 blinks correct!
    [thisexp.ftimes, thisexp.ftimesoff] = loadFrametimes(ftpath, thisexp.stimPara, expinfo.stimulusnames{ii}, ...
        expinfo.screen.delay, ksrastflag);
    
    
    disp('loading rasters');
    if ksrastflag   % update since using kilosort
        thisexp = loadKiloSortRasters( thisexp, ii, ksdata);
        
    elseif igorrasflag
        try
            thisexp = loadRasters( thisexp, ii);
        catch ME
            disp(ME.message);
            warning(['no spikes or ftimes found for experiment',num2str(ii),', check your experiment data!']);
        end
    end
    thisexp.info = expinfo;
    

    if ksrastflag
        thisexp.info.array.channelmap = ksdata.sort_params.channel_map;
        thisexp.info.array.channelpositions = ksdata.sort_params.channel_positions;
        thisexp.info.spikesorting.params_py = ksdata.sort_params.params_py;
        thisexp.info.spikesorting.stim_start_end = ksdata.sort_params.stim_start_end;

        if meaflag % to get MEA coordinates
            thisexp.info.array.meacoordinates = meacoord;
            m = meacoord;
            m(reshape(~ismember(1:size(m,1)*size(m,2),thisexp.mea.channelmap),size(m,1),size(m,2))) = 0;
            thisexp.info.array.meacoordmap = m;
        end
    else
        if meaflag,        thisexp.info.array.meacoordinates = meacoord;    end
    end
    disp('saving data');
    parsave(datapath,expinfo.stimulusnames{ii},thisexp.date,thisexp);
    %save([dataPath,'/Data Analysis/Raw Data/',files{ii},' for Experiment on ',exp.date,'.mat'],'-struct','exp');
    disp(['pre-analyzing stimulus ',expinfo.stimulusnames{ii}, ' is finito']);
end

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function experiment = parsave(dp,filename,expdate,expinput)
experiment = expinput;
save([dp,'/Data Analysis/Raw Data/',filename,' for Experiment on ',expdate,'.mat'],'-struct','experiment');
end
