function [ newsten ] = projOut(sten, projdir)
%PROJOUT Returns the STEN with the direction dir projected out
%   Detailed explanation goes here

%newsten=sten;
%sten=sten.STEN;
projdir = projdir(:);
inProds = sten*projdir;

PROJ = zeros(size(sten)); %initialize projection matrix
    
for i=1:size(sten,1)
    
    PROJ(i,:) = inProds(i)*projdir/sqrt(sum(projdir.^2));
end

%newSTEN=sten-PROJ; %remove the projection from the stEn

newsten=sten - PROJ; %remove the projection from the stEn

%newsten.STEN=newSTEN;

end

