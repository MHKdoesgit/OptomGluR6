

function []=concatenateRasters(datapath,chnum)
%--------------------------------------------------------------------------
%select the text files
if nargin < 1
    fstr={'*.txt', 'Text files'};
else
    if nargin<2, chnum = []; end
    if isnumeric(chnum), chnum = num2str(chnum); end
    
    fstr={[datapath,filesep,'rasters',filesep,'*_SP_C',chnum,'*.txt']};
end
[filenames, pathname]=uigetfile(fstr,'Select the rasters you want joined','MultiSelect', 'on');

if isequal(filenames,0); error('No rasters selected! Program is terminating...');
else
    disp('Rasters selected, proceeding to loading...');
end
% first get the indices of last cluster +1 and all stimuli number to set
% the defualt values.
spknames = dir([pathname,filesep,'*.txt']);   spknames = {spknames.name};
expnum = max(str2double(regexp(spknames,'\d*','Match','once')));
clusnum = cell2mat(regexp(filenames,'.txt','start'))-2;
clus = zeros(size(filenames));
for ii = 1:length(filenames)
    clus(ii) = str2double(filenames{ii}(clusnum(ii):clusnum(ii)+1));
end
clusnum = max(clus)+1;
%--------------------------------------------------------------------------
prompt = {'Enter the target cluster:', 'Enter the number of stimuli:'};
answer = inputdlg(prompt,'Important values',1,{num2str(clusnum),num2str(expnum)});
clustNum=str2double(answer{1});
expN=str2double(answer{2});
%--------------------------------------------------------------------------
%read channel number from filename
cstr=filenames{1};
channelNum=cstr(strfind(cstr,'C')+1:strfind(cstr,'.')-3);

for expId = 1:expN
    
    allspikes=[];
    for cellId=1:numel(filenames)
        undinds=strfind(filenames{cellId},'_');
        rawname=filenames{cellId}(undinds(1):end);
        spikes = load([pathname num2str(expId) rawname] );
        allspikes=[allspikes;spikes];
    end
    
    newspikes = sort(allspikes);
    newfilename = sprintf( '%s%d_SP_C%s%02d.txt',...
        pathname, expId, channelNum, clustNum);
    dlmwrite (newfilename, newspikes' ,'precision','%.5f')
end

end
