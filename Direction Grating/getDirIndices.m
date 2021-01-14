

function [dsi, dsi_angle, dsi_pval, osi, osi_angle, osi_pval, dsiDist, osiDist] = getDirIndices(...
    countMatrix, allAngles, ntrials)
%GETDIRINDICES Summary of this function goes here
%   Input:
%       countMatrix: a nTrials x nAngles matrix, populated by direction data
%       allAngles: a nAngles x 1 array of all stimulus angles in radians
% This function was written by Dimos and I adapted it to my use on
% 06.02.2018.

dsiAngles=allAngles;
osiAngles=2*allAngles;

dsiDist=NaN(1, ntrials);
osiDist=NaN(1, ntrials);

[m,n]=size(countMatrix);
for i=1:ntrials
    
    SMx = reshape(randperm(m*n), m, n);
    countMatrixShuffled=countMatrix(SMx);
    angleRatesShuffled=mean(countMatrixShuffled);
    
    %get shuffled dsi
    vectorsumTrial=sum(exp(0+1i*dsiAngles).*angleRatesShuffled);
    dsiTrial = abs(vectorsumTrial)/sum(abs(angleRatesShuffled));
    dsiDist(i)=dsiTrial;
    
    %get shuffled osi
    vectorsumTrial2=sum(exp(0+1i*osiAngles).*angleRatesShuffled);
    osiTrial= abs(vectorsumTrial2)/sum(abs(angleRatesShuffled));
    osiDist(i)=osiTrial;
    
end

angleRates=mean(countMatrix);

vectorsum=sum(exp(0+1i*dsiAngles).*angleRates);
dsi = abs(vectorsum)/sum(abs(angleRates));  %abs to avoid negative cases after pfr subtraction
dsi_angle = angle(vectorsum);
dsi_pval=1-sum(dsi>sort(dsiDist))/ntrials;
if dsi_pval==0
    dsi_pval=1/ntrials;
end

vectorsum2=sum(exp(0+1i*osiAngles).*angleRates);
osi= abs(vectorsum2)/sum(abs(angleRates));
osi_angle=angle(vectorsum2)/2 + pi/2;

osi_pval=1-sum(osi>sort(osiDist))/ntrials;
if osi_pval==0
    osi_pval=1/ntrials;
end

end

