

function varargout = plotReceptfield(gauss, img,varargin)
%
%%% plotReceptfield %%%
%
%
% This function plots the receptive filed and the Gaussian fit of it.
% it is the plotting section of receptfield function from Christian.
%
% ===============================Inputs====================================
%
%   gauss : parameters of Gaussian fit.
%   img : original image data.
%
%================================Output====================================
%
%   plot : plot of original data, Gaussian fit and corresponding residum.
%
% written by Mohammad, 18.04.2016


[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
f = @(s) s.coeff*normalPDF(X,Y,s.mu,s.sigma)+s.shift;
smoothimg = conv2(img,f(gauss),'same');

figure('position',[500 150 1000 800]);
subplot(2,2,1)
surf(Y,X,img,'edgecolor','none');
title('Original STA'); 
axis tight;

subplot(2,2,2)
surf(Y,X,f(gauss),'FaceColor','interp');
title('Optimized Gaussian fit');   
axis tight;

subplot(2,2,3)
surf(Y,X,smoothimg,'edgecolor','interp');
title('Smoothed STA with Gaussian fit');    
axis tight;

subplot(2,2,4)
surf(Y,X,f(gauss)-img,'FaceColor','interp');
title('Corresponding residuum');    
axis tight;


end




