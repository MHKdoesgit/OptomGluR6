
function [ clusters ] = loadClusters( datapath, bestScore, worstScore, varargin)
%
%%% loadClusters %%%
%
%
% This function will extract good channels of the excel file from spike
% sorting. It first check if a mat file from excel file was created or not,
% if not it create a .mat file with same name as original excel file
% within the same folder.
% if yes, it just load the mat file as clusters variable.
%
%================================Inputs====================================
%
%   datapath: folder path of each experiment.
%
%================================Output====================================
%
%   clusters : List of all good cells with their good clusters.
%   LisofGoodCells.mat: mat file of relevant columns in the excel file
%                       saved in the datapth folder.
%
% inspired by same function from Fernando,
% Written by Mohammad, 25.07.2014
% update for multiple excel selection on 01.08.2016.
% added best and worst score options on 2.1.2017.
% added NaN fill option for other excel file types on 15.02.2017.
% added new option to directly read the ksrasters file without a need to
% make new excel file, this is done for kilosort output on 19.02.2020.

if nargin > 3, loadexcel = varargin{1}; else, loadexcel = false; end

if ~(exist([datapath,'\*.xlsx'],'file')) && exist([datapath,'/ksrasters'],'dir')
    % this is to create same naming format as the old excel file style
    mainfoldname = lower(extractBetween(datapath,'D:\','\201'));
    if isempty(mainfoldname)
        mainfoldname = lower(extractBetween(datapath,'D:\','\202'));
    end
    if ~isempty(mainfoldname)
        mainfoldname = [upper(mainfoldname{1}(3)),mainfoldname{1}(4:end)];
    end
    electnum = cell2mat(extractBetween(datapath,'_','MEA'));
    expname = strrep(['201',extractAfter(datapath,'201')],['_',electnum,'MEA\'],'');
    mfn = [mainfoldname,electnum,'electrode_sortmodelrating',strrep(expname,'_','')];
    mfn = strrep(mfn,'/','');       mfn = strrep(mfn,'\','');
    matfilename = [datapath,'\CellsList_',mfn,'.mat'];
    if iscell(matfilename), matfilename = cell2mat(matfilename); end
    
    if exist(matfilename,'file')
        ListofGoodCells = struct2array(load (matfilename));
        disp('loading good clusters');
    else
        ListofGoodCells = struct2array(load([datapath,'/ksrasters/ksrasters.mat'],'clusters'));
        sortingquality = (ListofGoodCells(:,3) == bestScore | ListofGoodCells(:,3) <= worstScore);
        % 3rd colums is the cluster quality
        ListofGoodCells = ListofGoodCells(sortingquality,:);
        % saving the excel file as mat for faster loading
        save(matfilename,'-v7.3','ListofGoodCells');
        disp('Nope!!, Reading ksrasters file');
    end
    
else
    % this part is for igor spike sorting or for having the excel file for
    % ratings
    if nargin == 0
        [exelfilename,datapath] = uigetfile('D:/*.xlsx','Select the Excel file, Pour Favor');
        if isempty (exelfilename),  warning ('myfc:inputchk','File Not Found!!,Ola,Ola??');    end
    end
    
    % get excel file
    exelfilename = dir([datapath,'\*.xlsx']);
    
    if length(exelfilename) > 1
        warning ('myfc:inputchk','Too many Excel files on the dance floor');
        sel = listdlg('PromptString','Select a Excel file:','SelectionMode','single','ListString',{exelfilename.name},'ListSize',[400 150]);
        exelfilename = exelfilename(sel).name;
    else
        if ~isempty(exelfilename)
            exelfilename = exelfilename.name;
        else
            exelfilename = strrep(datapath(regexp(datapath,'201'):end),'_','');
            exelfilename = [strrep(exelfilename,'\',''),'.xlsx'];
        end
    end
    
    % check if mat file exists
    matfilename = [datapath,'\CellsList_',exelfilename(1:end-5),'.mat'];%dir([datapath,'\CellsList_*.mat']);
    disp('checking the existance of .mat file');
    
    if exist(matfilename,'file') && not(loadexcel)
        ListofGoodCells = struct2array(load (matfilename));
        disp('loading good clusters');
    else
        exceldata = xlsread ([datapath,'\',exelfilename]);
        if sum(isnan(exceldata( (find(not(isnan(exceldata(:,1))),1,'first'):find(not(isnan(exceldata(:,1))),1,'last')),1))) > 0
            exceldata(:,1) = fill_nans(exceldata(:,1));
        end
        excelrange = (exceldata(:,6) == bestScore | exceldata(:,6) <= worstScore);
        ListofGoodCells = exceldata( excelrange,[1,5:size(exceldata,2)]);
        ListofGoodCells = ListofGoodCells(:,not(all(isnan(ListofGoodCells))));      % remove unnecessary columns
        % saving the excel file as mat for faster loading
        save(matfilename,'-v7.3','ListofGoodCells');
        disp('Nope!!, Reading exel file');
    end
end

clusters = ListofGoodCells;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function A = fill_nans(A)
% Replaces the nans in each column with
% previous non-nan values.

for ii = 1:size(A,2)
    I = A(1,ii);
    for jj = 2:size(A,1)
        if isnan(A(jj,ii))
            A(jj,ii) = I;
        else
            I  = A(jj,ii);
        end
    end
end
end
