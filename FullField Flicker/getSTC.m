

function [ covMat ] = getSTC(stimEnsemble, binnedSpikes, varargin)
%GETSTC Get the STC matrix of a stimulus ensemble
%   Input:
%       stimEnsemble:
%       binnedSpikes:
%   Output:
%       covMat:

if nargin<3, isCentered=1; else, isCentered=varargin{1}; end
if size(stimEnsemble,1) ~= size(binnedSpikes,1), stimEnsemble = stimEnsemble'; end

covMat = zeros(size(stimEnsemble,2),size(stimEnsemble,2),size(binnedSpikes,2));

for ii = 1:size(binnedSpikes,2)
    %Make sure that the ensemble is empty of stimuli with zero spikes
    stimEn=stimEnsemble(binnedSpikes(:,ii) > 0,:);
    spkbin=binnedSpikes(binnedSpikes(:,ii) > 0,ii);
    
    %remove STA if needed
    sta=spkbin'*stimEn/sum(spkbin);
    stimEn=bsxfun(@minus, stimEn, sta*isCentered);
    
    %calculate covariance matrix fast
    covMat(:,:,ii) = stimEn'*bsxfun(@times,stimEn, spkbin)/(sum(spkbin)-1);
end
if size(covMat,3) == 1, covMat = squeeze(covMat(:,:,1)); end


end


% function [ STC ] = getSTC(sten, param)
%
%     stEn=sten.STEN;
%     binSpikes=sten.binSpikes;
%     sta=getSTA(sten);
%     STC0=0;
%     STC1=0;
%
%     for i=1:size(stEn,1)
%         STC0=STC0+(stEn(i,:)-sta)'*(stEn(i,:)-sta)*binSpikes(i);
%         STC1=STC1+(stEn(i,:))'*(stEn(i,:))*binSpikes(i);
%     end
%     STC0=STC0/(sum(binSpikes)-1);
%     STC1=STC1/(sum(binSpikes)-1);
%
%     if param==0
%         STC=STC0;
%     else
%         STC=STC1;
%     end
%
% end

