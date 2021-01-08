

function [c] = getEllipseFromParams(params,nsigma,numpoints)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<3; numpoints=50; end


gauss = rfroutines.getGaussFromParams(params);
c     = rfroutines.getEllipse(gauss, nsigma, numpoints); 
end

