

function plotchirpOptoLAP4(dp)

chipfolders = dir([dp,filesep,'Data Analysis',filesep,'*Chirp_Analysis*']);

savingpath = [dp,'/Data Analysis/Response_Comparison/Chirp_Analysis/'];
if not(exist(savingpath,'dir')), mkdir(savingpath); end

dat = cell(size(chipfolders,1),1);
for ii =1:size(chipfolders,1)
    chp = [chipfolders(ii).folder,filesep,chipfolders(ii).name];
    chpdat = dir([chp,filesep,'*Chirpstimulus_analysis_for_experiment_on*']);
    dat{ii} = load([chp,filesep,chpdat.name]);
end

dat = cell2mat(dat);

nstims = size(dat,1);
clus = dat(1).clusters;

ncolumns = 2;
nrows = nstims;
cols = lines(20);
cols = cols(6:end,:);
tt = 0.035;
titrset = {'Mesopic','Mesopic + L-AP4','Bleaching + L-AP4','Bleaching + Wash-out'};
msg = [];

for ii = 1:size(clus,1)
    
    h = figure('pos',[20 50 1750 850],'color','w','vis','off');
    for jj = 1:nstims
        
        [x,y] = raslin(dat(jj).rasters);
        p = dat(jj).para;
        
        subplot_tight(nrows,ncolumns,2*jj-1,tt)
        line(x{ii,end},y{ii,end},'color',cols(jj,:),'linewidth',0.25);
        if jj ==1
            hold on;
           % imagesc(linspace(0,p.changepoints(end),length(p.chirpstim)),p.Nrepeats+2:p.Nrepeats+2,repmat(p.chirpstim,2,1));
           % colormap(gray);
            plot(linspace(0,p.changepoints(end),length(p.chirpstim)),p.chirpstim + p.Nrepeats+2,'color','k');
            axis([0 round(p.changepoints(end)), 0 p.Nrepeats+5]);
        else
           axis([0 round(p.changepoints(end)), 0 p.Nrepeats+1]);  
        end
        
        ax = gca;          yticks(0:(p.Nrepeats+1)/2:p.Nrepeats+1);
        if jj~=4, ax.XColor = 'none';  end
        ax.TickLength= [0.005 0.005];           pbaspect([8 1 1]);      ylabel('trials');    xlabel('time (s)');
        pbaspect([5 1 1]);      title(titrset{jj});
        
        subplot_tight(nrows,ncolumns,2*jj,tt)
        
        line(dat(jj).psthtimes.psth,dat(jj).psth(ii,:),'color',cols(jj,:));
        xlim([0 round(p.changepoints(end))]);
        xticks(round(p.changepoints));      box off;
        ax = gca;       ax.TickLength= [0.005 0.005];  if jj~=4, ax.XColor = 'none';  end
        line([p.changepoints(2:end);p.changepoints(2:end)],[0;ax.YLim(2)],'color','r','linestyle','--');
%         if jj ==1
%         line(dat(jj).psthtimes.psth,p.chirpstim(1:3200)*(0.05*ax.YLim(2)) + ax.YLim(2)+5,'color','k');
%         end
        pbaspect([5 1 1]);              xlabel('time (s)');         ylabel('rate');
        title([titrset{jj}, '  (Rsq overall: ',num2str(round(dat(jj).rsqall(ii,end),2)),')'],'FontSize',9)
        
        
        
    end
    
    if isfield(p,'sortinfo'), chinfo =  p.sortinfo(ii); else, chinfo = clus(ii,:); end
    [filename,pngfilename] = rgcname('Chirp Stimulus', chinfo, p.date, ii);
    
    suptitle(h, filename,1);
    
    savepngFast(h, savingpath, pngfilename);
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', ii,size(clus,1),toc);
    fprintf(msg);
    close(h);

    
    
end



end


function [x,y] = raslin(r)

spkheight = 0.45;
[x,y] = deal(cell(size(r)));
for ii = 1: size(r,1)
    for jj = 1: size(r,2)
        rasresp = r{ii,jj};     %CelltoMatUE(r(jj,:));
        spkloc = repmat(1:size(rasresp,1),1,size(rasresp,2));
        spkoneline = [rasresp(:),rasresp(:),nan(numel(rasresp(:)),1)]';
        tboneline = [spkloc(:)-spkheight,spkloc(:)+spkheight, nan(numel(spkloc(:)),1)]' ; % + rasloc(jj)
        
        x{ii,jj} = spkoneline(:);
        y{ii,jj} = tboneline(:);
    end
end

end
