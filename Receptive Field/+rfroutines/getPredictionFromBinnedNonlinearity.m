

function [ preds ] = getPredictionFromBinnedNonlinearity(gens, centers, values)
%GETPREDICTIONFROMBINNEDNONLINEARITY Summary of this function goes here
%   Detailed explanation goes here

if any(isnan(values))
    preds = NaN(size(gens));
    return;
end

if any(diff(centers)==0) || any(diff(values)==0) % to avoid issues with noisy measurements
   preds = NaN(size(gens));
    return;
end

try
preds =  interp1(centers,values,gens,'linear','extrap');
catch ME
    disp(ME.message)
    preds =  interp1(sort(centers),sort(values),gens,'linear','extrap');
end


end

