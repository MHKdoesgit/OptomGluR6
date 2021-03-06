

function [ sigVecs, sigEigs, finalInts, finalMin, finalMax] = findSigEigs(STEnsemble, binnedspks, repTimes, confidence, n, varargin)
%FINDSIGEIGS Returns the significant eigenvalues and associated eigenvectors
%   n is the break parameter and sets an upper limit to the eigenvalue
%   detection.
%--------------------------------------------------------------------------
% fprintf('Checking significant eigenvalues and associated eigenvectors... '); tic;

sigVecs=[];
sigEigs=[];

if nargin > 5
    param = varargin{1};
else
    param = 1;
end

if n==0
    n=size(STEnsemble,2);
end

if isrow(binnedspks), binnedspks = binnedspks'; end

nStimuli=size(STEnsemble,1);

relDims=size(STEnsemble,2); %initialize the relevant dimensions
relDimsPos=1:relDims;

%trueSTC=getSTC(truesten, param); %get the true STC matrix
trueSTC= getSTC(STEnsemble,binnedspks, param); %get the true STC matrix
[trueVecs,trueEigs] = eig(trueSTC); %get the true eigenvectors-eigenvalues
trueEigs = sort(diag(trueEigs));

dumsten= STEnsemble;
dumSTEN = STEnsemble;

projlengths=zeros(nStimuli,1);

for i=1:nStimuli
    
    % proj=0;
    baseVecs = trueVecs(:,relDimsPos);
    inProds = dumSTEN(i,:)*baseVecs;
    
    %     for j=1:relDims
    %         proj=proj + inProds(j) * baseVecs(:,j);
    %     end
    %
    proj = sum(inProds .* baseVecs, 2);
    projlengths(i) = sqrt(sum(proj.^2));
end
shuffleEigs = zeros(repTimes, relDims,'gpuArray');

for i= 1 : repTimes
    
    shiftSTEN = randn(nStimuli, relDims,'gpuArray');
    %     for j=1:nStimuli
    %         shiftSTEN(j,:)=projlengths(j)*shiftSTEN(j,:)/sqrt(sum(shiftSTEN(j,:).^2));
    %     end
    
    shiftSTEN = bsxfun(@times, shiftSTEN./sqrt(sum(shiftSTEN.^2,2)), projlengths);
    
    shiftSTC = getSTCgpu(shiftSTEN,binnedspks);
    shiftEigs = sort(eig(shiftSTC));
    shuffleEigs(i,:) = shiftEigs';
end

trueInts(1, :) = quantile(gather(shuffleEigs),confidence/2);
trueInts(2, :) = quantile(gather(shuffleEigs),1-confidence/2);

%initialize indexes of stuff to check
iMin = find(trueEigs == min(trueEigs));
iMax = find(trueEigs == max(trueEigs));

dumInts = trueInts;
dumIntMin = min(trueInts(1,:));
dumIntMax = max(trueInts(2,:));

% nEigs=relDims-1; %number of eigenvalues to check from each side
% minEigs=trueEigs(iMin:nEigs);
% maxEigs=trueEigs(iMax-nEigs:iMax);

minEigs = trueEigs;
maxEigs = trueEigs;

counter=1;
% gfun1 = growdata2('rows');
% gfun2 = growdata2('columns');

while minEigs(1)<dumIntMin || maxEigs(end)>dumIntMax
    
    mmInd=find([dumIntMin-minEigs(1) maxEigs(end)-dumIntMax]==max([dumIntMin-minEigs(1) maxEigs(end)-dumIntMax]));
    
    if mmInd==1
        ind=find(trueEigs==minEigs(1));
        minEigs(1)=[];
        iMin=iMin+1;
        relDimsPos(1)=[];
    else
        ind=find(trueEigs==maxEigs(end));
        maxEigs(end)=[];
        iMax=iMax-1;
        relDimsPos(end)=[];
    end
    
    bigEig=trueEigs(ind);
    bigVec=trueVecs(:,ind);
    
    sigEigs=[sigEigs ; bigEig];     %#ok
    sigVecs=[sigVecs bigVec];       %#ok
    %     gfun1(bigEig);
    %     gfun2(bigVec);
    
    dumsten=projOut(dumsten, bigVec); %remove the projection from the stEn
    
    relDims=relDims-1;
    
    projlengths=zeros(nStimuli,1);
    for i=1:nStimuli
        %proj=0;
        baseVecs=trueVecs(:,relDimsPos);
        inProds=dumSTEN(i,:)*baseVecs;
        
        %         for j=1:relDims
        %             proj=proj+inProds(j)*baseVecs(:,j);
        %         end
        proj = sum(inProds .* baseVecs, 2);
        projlengths(i)=sqrt(sum(proj.^2));
    end
    
    % shiftsten=dumsten;
    shuffleEigs=zeros(repTimes, relDims, 'gpuArray');
    
    for i=1:repTimes
        shiftSTEN=randn(nStimuli, relDims, 'gpuArray');
        %         for j=1:nStimuli
        %             shiftSTEN(j,:)=projlengths(j)*shiftSTEN(j,:)/sqrt(sum(shiftSTEN(j,:).^2));
        %         end
        shiftSTEN = bsxfun(@times, shiftSTEN./sqrt(sum(shiftSTEN.^2,2)), projlengths);
        shiftSTC = getSTCgpu(shiftSTEN,binnedspks, param);
        shiftEigs=sort(eig(shiftSTC));
        shuffleEigs(i,:)=shiftEigs';
    end
    
    dumInts= [quantile(gather(shuffleEigs),confidence/2); quantile(gather(shuffleEigs),1-confidence/2)];
    dumIntMin=min(dumInts(1,:));
    dumIntMax=max(dumInts(2,:));
    
    if counter==n
        break ;
    end
    counter=counter+1;
end

% sigEigs=gfun1();
% sigVecs=gfun2();

finalMin=iMin;
finalMax=iMax;
finalInts=dumInts;

% fprintf('Done! Took %2.2f s...\n', toc);

end



function [ covMat ] = getSTCgpu(stimEnsemble, binnedSpikes, varargin)

if size(stimEnsemble,1) ~= size(binnedSpikes,1), stimEnsemble = stimEnsemble'; end
covMat = zeros(size(stimEnsemble,2),size(stimEnsemble,2),size(binnedSpikes,2),'gpuArray');

for ii = 1:size(binnedSpikes,2)
    %Make sure that the ensemble is empty of stimuli with zero spikes
    stimEn=stimEnsemble(binnedSpikes(:,ii) > 0,:);
    spkbin=binnedSpikes(binnedSpikes(:,ii) > 0,ii);
    %remove STA if needed
    sta     =   spkbin'*stimEn/sum(spkbin);
    stimEn  =   bsxfun(@minus, stimEn, sta);
    %calculate covariance matrix fast
    covMat(:,:,ii) = stimEn'*bsxfun(@times,stimEn, spkbin)/(sum(spkbin)-1);
end
if size(covMat,3) == 1, covMat = squeeze(covMat(:,:,1)); end

end

