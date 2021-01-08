

function [ experiment ] = loadExperimentData( datapath , varargin)
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
%    dataPath : path to experiment folder
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

if nargin < 1
   dataPath = uigetdir();
end


if (any(datapath == 0) || isempty(datapath) || ~isfolder(datapath))
    error('You must supply a dataPath path');
end

bestScore = 1;  % by default only the cells with score 1 are considerd as good.
worstScore = 1;

if nargin == 2  % if only one input is given that is considered as worst score and 1 as best score.
    bestScore = 1;
    worstScore = varargin{1};
end

if nargin == 3    % both best and worst scores are defined by user.
    bestScore = varargin{1};
    worstScore = varargin{2};
end

%lightProjection = questdlg('Choose your projector','Projection System','LightCrafter','OLED','LightCrafter');
[recordingSetup, lightProjection] = getrecordingsetupname();
switch lower(lightProjection)
    case 'lightcrafter'
        screenWidth = 864;
        screenHeight = 480;
        monitorDelay = 35;  % 36ms for uv and 34 ms for green
    case 'oled'
        screenWidth = 800;
        screenHeight = 600;
        monitorDelay = 25;
    otherwise
        error ('No light projection device is selected, please choose one!');
end

% get screen refresh rate
q = questdlg('select screen refresh rate (Hz):','refresh rate','60 Hz','75 Hz','85 Hz','60 Hz');
screenrefreshrate = str2double(extractBefore(q,' Hz'));

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

experiment.originalFolder = datapath;
[~, expName] = fileparts(experiment.originalFolder);
experiment.expName = expName;
experiment.date = datemaker(datapath);

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
if not(exist(paramFile, 'file'))
    parameterstxtfromftimes(ftpath(1:end-1));
end

disp('loading experiment parameters');
[files, num, header] = loadParameters(paramFile);
experiment.parameters.header = header;
experiment.parameters.files = files;
experiment.parameters.numFiles = num;

disp('loading good cells');
experiment.clusters = loadClusters( datapath ,bestScore, worstScore);

disp('get MEA coordinates');
try
    meacoord = getMEAcoordinatesfromMCD(datapath);
    meaflag= true;
catch ME
    disp(ME.message)
    warning('There aint no MEA pictures and MEA coordinates files here!');
    meaflag = false;
end

disp('loading stimuli parameters');

for ii = 1:numel(files)
    
    thisexp = experiment;
    thisexp.stimPara = loadStimulusParameters(files{ii}, screenWidth, screenHeight);
    
    disp('loading frametimes'); % can be done cleaner but necessary to have 1 blinks correct!
    if ksrastflag && any(isfield(thisexp.stimPara,{'nblinks','nblink','Nblink','Nblinks',}))
        stimuliFrames = load([ftpath, files{ii}, '_frametimings.mat']);
        if mod(thisexp.stimPara.nblinks,2) == 1 % in case of nblink 1 or odd numbers
            stimuliFrames.ftimes = sort([stimuliFrames.ftimes;stimuliFrames.ftimesoff],'ascend')';
        end
    else
        stimuliFrames = load([ftpath, files{ii}, '_frametimings.mat'], 'ftimes','ftimesoff');
    end
    % 25-35 ms monitor delay is added to frametimes
    if ksrastflag
        thisexp.ftimes = (stimuliFrames.ftimes + (monitorDelay/1e3)); % frametime in ksrasters in already in secs!
        thisexp.ftimesoff = (stimuliFrames.ftimesoff + (monitorDelay/1e3)); % frame offset
    else
        if isfield(stimuliFrames,'ftimesoff') % this for running igor with new frametimes
            if any(isfield(thisexp.stimPara,{'nblinks','nblink','Nblink','Nblinks',}))
                if mod(thisexp.stimPara.nblinks,2) == 1 % in case of nblink 1 or odd numbers
                    stimuliFrames.ftimes = sort([stimuliFrames.ftimes;stimuliFrames.ftimesoff],'ascend')';
                end
            end
            thisexp.ftimes = (stimuliFrames.ftimes + (monitorDelay/1e3)); % frametime in ksrasters in already in secs!
            thisexp.ftimesoff = (stimuliFrames.ftimesoff + (monitorDelay/1e3)); % frame offset
        else
            thisexp.ftimes = (stimuliFrames.ftimes + monitorDelay) / 1e3; % ftimes saved in ms, while spikes are in s!
        end
    end
    
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
    
    thisexp.screen          =    [screenWidth screenHeight];
    thisexp.lightprojection =    lightProjection;
    thisexp.refreshrate     =    screenrefreshrate;
    thisexp.recordingsetup  =    recordingSetup;
    thisexp.samplingrate    =    double(header(2));
    
    thisexp.type            =   'MEA';
    if meaflag
        thisexp.meacoordinates = meacoord;
    end
    
    if ksrastflag
        thisexp.mea.numeletrodes = header(1);
        thisexp.mea.channelmap = ksdata.sort_params.channel_map;
        thisexp.mea.channelpos = ksdata.sort_params.channel_positions;
        if meaflag % to get MEA coordinates
            thisexp.mea.meacoordinates = meacoord;
            m = meacoord;
            m(reshape(~ismember(1:size(m,1)*size(m,2),thisexp.mea.channelmap),size(m,1),size(m,2))) = 0;
            thisexp.mea.meacoordmap = m;
        end
    else
        if meaflag,        thisexp.meacoordinates = meacoord;    end
    end
    disp('saving data');
    parsave(datapath,files{ii},thisexp.date,thisexp);
    %save([dataPath,'/Data Analysis/Raw Data/',files{ii},' for Experiment on ',exp.date,'.mat'],'-struct','exp');
    disp(['pre-analyzing stimulus ',files{ii}, ' is finito']);
end

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function experiment = parsave(dp,filename,expdate,expinput)
experiment = expinput;
save([dp,'/Data Analysis/Raw Data/',filename,' for Experiment on ',expdate,'.mat'],'-struct','experiment');
end

%--------------------------------------------------------------------------------------------------%

function [recordingSetup, lightSource] = getrecordingsetupname()

setupnames = {'Aragorn-lightcrafter','Aragorn-OLED','Bilbo-lightcrafter','Bilbo-OLED','Elrond-OLED',...
    'Gandalf-lightcrafter','Gandalf-OLED','Gandalf-CMOS-lightcrafter','Gandalf-CMOS-OLED',...
    'Legolas-lightcrafter','Legolas-OLED','Bombadil-OLED','Saruman-OLED','Lightcrafter','OLED','Unknown'};

exprecsetup = listdlg('ListString',setupnames,'PromptString','Select the recording setup:',...
    'SelectionMode','single','ListSize',[180 300]);

recordingSetup = setupnames{exprecsetup};
if contains(recordingSetup,{'-lightcrafter','-OLED'})
    lightSource = extractAfter (recordingSetup,'-');
    lightSource = [upper(lightSource(1)),lightSource(2:end)];
    recordingSetup = extractBefore (recordingSetup,'-');
else
    lightSource = setupnames{exprecsetup};
    recordingSetup = 'Unknown';
end

end
