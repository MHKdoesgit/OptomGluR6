

function [ STA ] = getSTA(stimEnsemble, binnedSpikes)
%GETSTA Returns the non-normalized STA
%   Detailed explanation goes here

% STA=0;
% STEN=sten.STEN;
% binSpikes=sten.binSpikes;
% 
% for i=1:size(STEN,1)
%    STA=STA+STEN(i,:)*binSpikes(i);
% end
% 
% STA=STA/sum(binSpikes);


%GETSTA Returns the non-normalized STA
%   Input:
%       
%       binnedSpikes:
%   Output:
%       STA:

STA=binnedSpikes*stimEnsemble/sum(binnedSpikes);

end

