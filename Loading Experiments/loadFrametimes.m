

function [ftimes, ftimesoff] = loadFrametimes(ftpath, stimpara, stimname, monitordelay, ksrastflag, varargin)
%
%%% loadFrametimes %%%
%
%
% This function loads the frametimes for the mat file created by the spike
% sorting code. It also adds the monitor delay to each frametimes file and
% for the stimuli with 1 blinks, it joins the ftimes on and off together to
% make 1 blink ftimes.
%
% ===============================Inputs====================================
%
%    ftpath : path to frametimes folder
%    stimpara : stimulus parameters to get the nblinks.
%    stimname : name of the file for loading the correct frametimings.
%    monitordelay : monitor dealy in ms.
%    ksrastflag : to load new format of framtiming with ftimes on and off.
%
%================================Output====================================
%
%   ftimes : onset of each pulse in seconds.
%   ftimesoff : offset of each pulse in seconds.
%
% written by Mohammad on 07.01.2021.

if ksrastflag && any(isfield(stimpara,{'nblinks','nblink','Nblink','Nblinks'}))
    fn = {'nblinks','nblink','Nblink','Nblinks'};
    fn = fn{isfield(stimpara,{'nblinks','nblink','Nblink','Nblinks'})};
    stimuliFrames = load([ftpath, stimname, '_frametimings.mat']);
    if mod(stimpara.(fn),2) == 1 % in case of nblink 1 or odd numbers
        stimuliFrames.ftimes = sort([stimuliFrames.ftimes;stimuliFrames.ftimesoff],'ascend')';
    end
else
    stimuliFrames = load([ftpath, stimname, '_frametimings.mat'], 'ftimes','ftimesoff');
end
% 25-35 ms monitor delay is added to frametimes
if ksrastflag
    ftimes = (stimuliFrames.ftimes + (monitordelay/1e3)); % frametime in ksrasters in already in secs!
    ftimesoff = (stimuliFrames.ftimesoff + (monitordelay/1e3)); % frame offset
else
    if isfield(stimuliFrames,'ftimesoff') % this for running igor with new frametimes
        if any(isfield(stimpara,{'nblinks','nblink','Nblink','Nblinks'}))
            fn = {'nblinks','nblink','Nblink','Nblinks'};
            fn = fn{isfield(stimpara,{'nblinks','nblink','Nblink','Nblinks'})};
            if mod(stimpara.(fn),2) == 1 % in case of nblink 1 or odd numbers
                stimuliFrames.ftimes = sort([stimuliFrames.ftimes;stimuliFrames.ftimesoff],'ascend')';
            end
        end
        ftimes = (stimuliFrames.ftimes + (monitordelay/1e3)); % frametime in ksrasters in already in secs!
        ftimesoff = (stimuliFrames.ftimesoff + (monitordelay/1e3)); % frame offset
    else
        ftimes = (stimuliFrames.ftimes + monitordelay) / 1e3; % ftimes saved in ms, while spikes are in s!
    end
end

if isrow(ftimes), ftimes = transpose(ftimes); end
if isrow(ftimesoff), ftimesoff = transpose(ftimesoff); end

end