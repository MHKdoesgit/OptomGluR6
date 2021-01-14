

function [sten]=getSTEN(stimulus, spt, windowN)
%GETSTEN Computes the stimulus ensemble matrix given the spike train and the window of points. 
%The STEN is a struct with 2 elements: STEN is the spike
%triggered ensemble matrix (N x p), where N is the number of
%spike-triggered stimuli and p the correscponding time points and,
%binSpikes (N x 1) is the number of spikes that each stimulus evoked.
%I re-wrote this function to better suit the data structure. 
% written by Dimos, 8 March 2016.

stim=stimulus(:);
spt(1:windowN) = 0; % remove spikes that occur in 1st timeWindow time (no stimulus for whole block of time prior to those spikes)
sp_ind=find(spt);
ens=zeros(length(sp_ind), windowN);
n_spikes= zeros(length(sp_ind),1); % number of spikes per stimulus
idx=1;

for i=sp_ind
   ens(idx, :) = stim(i-windowN+1:i);
   n_spikes(idx)=spt(i);
   
   idx=idx+1;
end

sten=struct();
sten.STEN=ens;
sten.binSpikes=n_spikes;

end