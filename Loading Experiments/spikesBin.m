function [ spikesbin ] = spikesBin( ftimes, spikeTimes)
%
%%% spikesBinner %%%
%
%
% This function bin the spiketimes bsed on stimulus ftimes.
% for alternative method check spikesBinner_V2.
%
% ===============================Inputs====================================
%
%    ftimes : stimulus frame timings. 
%    spikeTimes : Experiment spike timings.
%
%================================Output====================================
%
%   spikesbin : binned spiketimes.
%   
% Note : Note this function cannot bin spikes at resolution below one frame
%        in that case use spikesBinner instead.
%
% Based on code by Michael Weick,
% Copied from spikeBin function from Fernando,
% modified by Mohammad, 30.07.2014



if (numel(ftimes) < 2)
    spikesbin = numel(spikeTimes);
    return
end


lastFrame = numel(ftimes); % each 2nd frame puts a new pulse
avgFrameInt = mean(diff(ftimes));


% Set t = 0 to be the first frame and remove the spikes which didn't
% happen during stimulus presentation.
useSpikes = spikeTimes - ftimes(1);
useSpikes(useSpikes <= 0) = []; % before stimulus
useSpikes(useSpikes > (lastFrame * avgFrameInt )) = [];

%useSpikes(useSpikes > ((lastFrame + depth)*avgFrameInt) ) = [];

binsPerFrame =round( 1000 * avgFrameInt );
timePerBin = avgFrameInt/binsPerFrame;

spksBin = zeros( lastFrame * binsPerFrame, 1);
spksBin( ceil(useSpikes / timePerBin) ) = 1;

% Nice trick to "rebin" spikes per frame:
spikesbin = sum(reshape( spksBin, [ binsPerFrame lastFrame ] ), 1);


end

