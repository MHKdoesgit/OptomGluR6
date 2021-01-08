

function stimnewnames = zerobeforeStimNames(foldernamepath)
%
%%% zerobeforeStimNames %%%
%
%
% This function add one zero before the numbers 1 to 9 of any input name to
% have them sorted easier in Windows. for example 1_st becomes 01_st
%
%================================Inputs====================================
%
%   foldernamepath : folder path.
%
%================================Output====================================
%
%   stimnewnames : new name list and as they were renamed in the folder.
%
% written by Mohammad, 21.04.2020.

tic;
stimoldnames = dir(foldernamepath);
stimoldnames = {stimoldnames(3:end).name};

stimnums = extractBefore(stimoldnames,'_');

stimnumsdig = cellfun(@str2double,stimnums);

num2change = find(stimnumsdig < 10);

stimnewnames = stimoldnames;
for ii = 1:length(num2change)
    
    stimnewnames{num2change(ii)} = [num2str(stimnumsdig(num2change(ii)),'%02d'),'_',extractAfter(stimoldnames{num2change(ii)},'_')];
end


for jj = 1 : length(stimoldnames)
    if not(strcmp(stimoldnames{jj}, stimnewnames{jj}))
        movefile( fullfile(foldernamepath, stimoldnames{jj}), fullfile(foldernamepath, stimnewnames{jj}) );
    end
end
toc;

end