

function [ nly, nlx, nlystd] = getDichromaticNonlinearities(NLtype, linOutputCol1, linOutputCol2, spikes, nbins, dt, varargin)

if nargin > 6, conditionIdx = varargin{1}; else, conditionIdx = 0.1; end

if size(spikes,1)~= size(linOutputCol1,1), spikes = spikes'; end % it has to be colums to work with accumarray
Ncells = size(linOutputCol1,1);

switch lower(NLtype)
    
    case {'overall','bothcolors','together','combined','all','o'}
        
        [nly, nlx, nlystd] = deal(zeros(Ncells, nbins,'single'));
        for ii = 1: Ncells
            celllinoutcol1 = linOutputCol1(ii,:);
            celllinoutcol2 = linOutputCol2(ii,:);
            if all(isnan(sum([celllinoutcol1(:),celllinoutcol2(:)]))); continue; end
            [ nly(ii,:), nlx(ii,:), nlystd(ii,:)] = overallNL(celllinoutcol1, celllinoutcol2, spikes(ii,:), nbins, dt);
        end
        
    case {'conditional','cond','c'}
        
        [nly, nlx, nlystd] = deal(zeros(Ncells, 2, nbins,'single'));
        for ii = 1: Ncells
            celllinoutcol1 = linOutputCol1(ii,:);
            celllinoutcol2 = linOutputCol2(ii,:);
            if all(isnan(sum([celllinoutcol1(:),celllinoutcol2(:)]))); continue; end
            [ nly(ii,:,:), nlx(ii,:,:), nlystd(ii,:,:)] = conditionalNL(celllinoutcol1, celllinoutcol2, spikes(ii,:), nbins, dt , conditionIdx);
        end
        
    case {'2d','nl2d', '2dimentional','2-d','2-dimentional','twodim'}
        
        [nly, nlystd] = deal(zeros(Ncells, nbins, nbins,'single'));
        nlx = zeros(Ncells, 2, nbins,'single');
        for ii = 1: Ncells
            celllinoutcol1 = linOutputCol1(ii,:);
            celllinoutcol2 = linOutputCol2(ii,:);
            if all(isnan(sum([celllinoutcol1(:),celllinoutcol2(:)]))); continue; end
            [ nly(ii,:,:), nlx(ii,:,:), nlystd(ii,:,:)] = twoDimNL(celllinoutcol1, celllinoutcol2, spikes(ii,:), nbins, dt);
        end
        
        % for fun plotting
        %for ii = 1:Ncells
        %    surf(squeeze(nlx(ii,1,:)),squeeze(nlx(ii,2,:)),...
        %        smoothdata(squeeze(nly(ii,:,:)),'gaussian',5),'EdgeColor','texturemap','FaceColor','interp')
        %    pause;
        %end
        
    case {'marginal','marg','m'}
        
        [nly, nlx, nlystd] = deal(zeros(Ncells, 2, nbins,'single'));
        for ii = 1: Ncells
            celllinoutcol1 = linOutputCol1(ii,:);
            celllinoutcol2 = linOutputCol2(ii,:);
            if all(isnan(sum([celllinoutcol1(:),celllinoutcol2(:)]))); continue; end
            [ nly(ii,:,:), nlx(ii,:,:), nlystd(ii,:,:)] = marginalNL(celllinoutcol1, celllinoutcol2, spikes(ii,:), nbins, dt );
        end
        
    otherwise
        error('the fuck is this nonlinearity type! choose a proper one!');
end

end

function [ nly, nlx, nlystd] = overallNL(gensigCol1, gensigCol2, spk, nbins, dt)

qtStep = 1/nbins;
qtEdges = 0:qtStep:1;
% we combine the two colors to get overall nonlinearities
linearOutput = gensigCol1 + gensigCol2;

genEdges = quantile(linearOutput, qtEdges);

[bincnts, ~, spikebins] = histcounts(linearOutput,genEdges);

nly   = accumarray(spikebins', spk, [nbins 1], @sum);
sdvalues = accumarray(spikebins', spk, [nbins 1], @std);

nly   = nly./bincnts'/dt;
nlystd = sdvalues./sqrt(bincnts')/dt;

nlx = genEdges(1:nbins)+diff(genEdges)/2;

end


function [ nly, nlx, nlystd] = conditionalNL(gensigCol1, gensigCol2, spk, nbins, dt , conditionIdx)

qtStep = 1/nbins;
qtEdges = 0:qtStep:1;
% do conditioning
idxcol1 = (abs(gensigCol2) <= conditionIdx); % index in the other color
idxcol2 = (abs(gensigCol1) <= conditionIdx);

genEdgesCol1 = quantile(gensigCol1(idxcol1), qtEdges);
genEdgesCol2 = quantile(gensigCol2(idxcol2), qtEdges);

if all(isnan(genEdgesCol1)) || all(isnan(genEdgesCol2))
    [nly, nlx, nlystd] = deal(nan(2, nbins));
    fprintf('no response in the conditioning range, color 1: %.4f, color  2: %.4f\n',...
        min(abs(gensigCol2)),min(abs(gensigCol1)));
    return;
end

[acol1, ~, spikebinsCol1] = histcounts(gensigCol1, genEdgesCol1);
[acol2, ~, spikebinsCol2] = histcounts(gensigCol2, genEdgesCol2);

% to avoid zero index with accumarray
[~,~,spikebinsCol1] = unique(spikebinsCol1);
[~,~,spikebinsCol2] = unique(spikebinsCol2);

valuesCol1   = accumarray(spikebinsCol1, spk, [nbins+1 1], @sum);
sdvaluesCol1 = accumarray(spikebinsCol1, spk, [nbins+1 1], @std);
valuesCol1   = valuesCol1(2:end);
sdvaluesCol1 = sdvaluesCol1(2:end);

valuesCol2   = accumarray(spikebinsCol2, spk, [nbins+1 1], @sum);
sdvaluesCol2 = accumarray(spikebinsCol2, spk, [nbins+1 1], @std);
valuesCol2   = valuesCol2(2:end);
sdvaluesCol2 = sdvaluesCol2(2:end);


valuesCol1   = valuesCol1./acol1'/dt;
sevaluesCol1 = sdvaluesCol1./sqrt(acol1')/dt;

valuesCol2   = valuesCol2./acol2'/dt;
sevaluesCol2 = sdvaluesCol2./sqrt(acol2')/dt;

centersCol1 = genEdgesCol1(1:nbins)+diff(genEdgesCol1)/2;
centersCol2 = genEdgesCol2(1:nbins)+diff(genEdgesCol2)/2;

% put outputs together, first row: first color, second row: second color
nly = [valuesCol1,valuesCol2]';
nlystd = [sevaluesCol1 , sevaluesCol2]';
nlx = [centersCol1 ; centersCol2];


end


function [ nly, nlx, nlystd] = twoDimNL(gensigCol1, gensigCol2, spk, nbins, dt)

qtStep = 1/nbins;
qtEdges = 0:qtStep:1;

genEdgesCol1 = quantile(gensigCol1, qtEdges);
genEdgesCol2 = quantile(gensigCol2, qtEdges);

[bincnts,~,~,spikebinsX,spikebinsY] = histcounts2(gensigCol1,gensigCol2,genEdgesCol1,genEdgesCol2);

values2d   = accumarray([spikebinsX;spikebinsY]', spk, [nbins nbins], @sum);
sdvalues2d   = accumarray([spikebinsX;spikebinsY]', spk, [nbins nbins], @std);

values2d   = values2d./bincnts/dt;
sevalues2d = sdvalues2d./sqrt(bincnts)/dt;

centersCol1 = genEdgesCol1(1:nbins)+diff(genEdgesCol1)/2;
centersCol2 = genEdgesCol2(1:nbins)+diff(genEdgesCol2)/2;

nly = values2d;
nlystd = sevalues2d;
nlx = [centersCol1 ; centersCol2];

end


function [ nly, nlx, nlystd] = marginalNL(gensigCol1, gensigCol2, spk, nbins, dt)

qtStep = 1/nbins;
qtEdges = 0:qtStep:1;

genEdgesCol1 = quantile(gensigCol1, qtEdges);
genEdgesCol2 = quantile(gensigCol2, qtEdges);

[acol1, ~, spikebinsCol1] = histcounts(gensigCol1, genEdgesCol1);
[acol2, ~, spikebinsCol2] = histcounts(gensigCol2, genEdgesCol2);

valuesCol1   = accumarray(spikebinsCol1', spk, [nbins 1], @sum);
sdvaluesCol1 = accumarray(spikebinsCol1', spk, [nbins 1], @std);

valuesCol2   = accumarray(spikebinsCol2', spk, [nbins 1], @sum);
sdvaluesCol2 = accumarray(spikebinsCol2', spk, [nbins 1], @std);

valuesCol1   = valuesCol1./acol1'/dt;
sevaluesCol1 = sdvaluesCol1./sqrt(acol1')/dt;

valuesCol2   = valuesCol2./acol2'/dt;
sevaluesCol2 = sdvaluesCol2./sqrt(acol2')/dt;

centersCol1 = genEdgesCol1(1:nbins)+diff(genEdgesCol1)/2;
centersCol2 = genEdgesCol2(1:nbins)+diff(genEdgesCol2)/2;

% put outputs together, first row: first color, second row: second color
nly = [valuesCol1, valuesCol2]';
nlystd = [sevaluesCol1 , sevaluesCol2]';
nlx = [centersCol1 ; centersCol2];

end
