function [ centGen, rateEst ] = probNlinEst(stimulus, sten, recFilters, nBins)
%PROBNLINEST Returns the nonlinearity estimates for each filter in recFilters
%   nBins is usually set (in papers) to 40

stEn=sten.STEN;
estGens=zeros(size(recFilters,1), length(stimulus));
qtStep=1/nBins;
qts=qtStep/2:qtStep:1-qtStep/2;

for i=1:size(recFilters,1)
    
    genS=conv(stimulus,fliplr(recFilters(i,:)), 'full');
    genS(length(stimulus)+1:end)=[];
    estGens(i,:)=genS;
    
end

rateEst=zeros(size(estGens,1),nBins);
centGen=zeros(size(estGens,1),nBins);

for i=1:size(estGens,1)
    
    estGen=estGens(i,:);
    centers = quantile(estGen,qts);
    
    h=hist(estGens(i,:), centers)/size(estGens,2);
    hs=hist(stEn*recFilters(i,:)', centers)/size(stEn,1);
    
    rateEst(i,:)=hs./h;
    centGen(i, :)=centers;
end


end

