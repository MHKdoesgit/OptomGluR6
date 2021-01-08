

function [ftnames, params] = parameterstxtfromftimes(datapath, numchannels, samplerate, paramsnumch)
%
%%% parameterstxtfromftimes %%%
%
%
% This function makes legacy parameters text files based on the names of
% the mcd files that is used in frametimes folder. This is an auxiliary
% function for new spike sorting with kilosort and used for backward
% compatibility with Igor spike sorting.
%
% ===============================Inputs====================================
%
%    datapath : path to frametimes folder created by KS2 code.
%    numchannels : number of electrodes to make the header file.(optional, gui).
%    samplerate : sampling rate used to make the header file.(optional, gui).
%    paramsnumch : this flag make an extra text file without the header
%    which was used to be used by old frametimes code. Not needed at all!
%
%================================Output====================================
%
%   ftnames : stimuli, frametime files and parameters names.
%   params : header of the parameters.txt, include numelectrodes, sampling
%            rates, etc.
%
% written by Mohammad, 19.02.2020 based for KS2 marmoset output.

if nargin < 2
    prompt = {'MEA electrode numbers:','Sampling rate:'};
    title = 'parameters from frametimes';
    dims = [1 50];
    definput = {'252','10000'};
    answer = inputdlg(prompt,title,dims,definput);
    numchannels = str2double(answer{1});
    samplerate = str2double(answer{2});
end

if nargin < 4, paramsnumch = false; end

ftfiles = dir([datapath,'/*_frametimings.mat']);
ftfiles = {ftfiles.name}';

ftnames = extractBefore(ftfiles,'_frametimings.mat');

% tab is char(9)
headerline = [num2str(numchannels),char(9),num2str(samplerate),char(9),'1',...
    char(9),num2str(size(ftnames,1)),char(9),'0',char(9),char(9),newline];

fid = fopen([fileparts(datapath),'/parameters.txt'],'W');
fwrite(fid,headerline);
for ii = 1:size(ftnames)
    %fwrite(fid, [ftnames{ii},char(9),char(9),'1',newline]);
    fprintf(fid,'%-65s%d\n', ftnames{ii}, 1);
end
fclose(fid);

if paramsnumch
    fid = fopen([fileparts(datapath),'/parameters',num2str(numchannels),'.txt'],'W');
    for ii = 1:size(ftnames)
        %fwrite(fid, [ftnames{ii},char(9),char(9),'1',newline]);
        fprintf(fid,'%-65s%d\n', ftnames{ii}, 1);
    end
    fclose(fid);
    
end

params = str2num(headerline); %#ok

end