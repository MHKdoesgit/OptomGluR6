

function CheckAmpBeforeAfterBleaching(dp, stimbefnum, stimafternum)


rawdatpath = [dp, filesep, 'Data Analysis',filesep,'Raw Data',filesep];

stimnames = dir([rawdatpath,'*_*.mat']);
stimnames = {stimnames.name}';

clus = struct2array(load([rawdatpath,stimnames{1}],'clusters'));
[allamps, allspks] = deal(cell(size(stimnames,1),size(clus,1)));
for ii = 1:size(stimnames,1)
    allamps(ii,1:size(clus,1)) = struct2array(load([rawdatpath,stimnames{ii}],'amplitudes'));
    allspks(ii,1:size(clus,1)) = struct2array(load([rawdatpath,stimnames{ii}],'spiketimes'));
end

stimbefore = cell2mat(stimnames(contains(stimnames,num2str(stimbefnum,'%02d_'))));
stimafter = cell2mat(stimnames (contains(stimnames,num2str(stimafternum,'%02d_'))));


stimbefdat = load([rawdatpath,stimbefore]);
stimaftdat = load([rawdatpath,stimafter]);

stimstartend = stimbefdat.info.spikesorting.stim_start_end / stimbefdat.info.samplingrate;
nstplt = 25000;
npall = 5000;
cols = lines(10);
tt = 0.035;

for ii = 1: size(clus,1)
%%

    ampbef = stimbefdat.amplitudes{ii};
    spkbef = stimbefdat.spiketimes{ii};
    
    ampaft = stimaftdat.amplitudes{ii};
    spkaft = stimaftdat.spiketimes{ii};
    
    maxdur = ceil(max([max(spkbef),max(spkaft)])/60);
    tbins = linspace(0,maxdur,maxdur-1);
    
    
    binnedspkbef = histcounts(spkbef/60,tbins);
    binnedspkaft = histcounts(spkaft/60,tbins);
    
    
if size(ampbef,1) > nstplt, ns = nstplt; else, ns = size(ampbef,1); end
   
befid = randperm(size(ampbef,1),ns);
aftid = randperm(size(ampaft,1),ns);

expdur =  ceil(max(stimstartend(:))/60);
expbins = linspace(0,expdur,expdur-1);
binnedallspks = histcounts(cell2mat(cellfun(@plus,allspks(:,ii),num2cell(stimstartend(:,1)),'un',0))/60, expbins);


h = figure('pos',[100 10 1500 1000],'color','w','vis','on');

subplot_tight(5,1,1,tt)
plot(spkbef(befid)/60, ampbef(befid),'.','color',cols(1,:));
yAx = ceil(max(ampbef(befid))/2)*2;
axis([0 maxdur 5 yAx]);
xticks(0:10:400);               yticks(5:5:50);         box off;
%xlabel('Time (min)');   
ylabel('Std from noise');
title(['Before bleaching, stimulus:',strrep(stimbefore(1:end-4),'_','-')]);
set(gca,'Xcolor','none');

subplot_tight(5,1,2,tt)
plot(spkaft(aftid)/60, ampaft(aftid),'.','color',cols(2,:));
axis([0 maxdur 5 yAx]);
xticks(0:10:400);               yticks(5:5:50);         box off;
xlabel('Time (min)');   ylabel('Std from noise');
title(['After bleaching, stimulus:',strrep(stimafter(1:end-4),'_','-')]);


subplot_tight(5,1,3,tt)

for jj = 1:size(allamps,1)
   if size(allamps{jj,ii},1) > npall, ns = npall; else, ns = size(allamps{jj,ii},1); end 
   id = randperm(size(allamps{jj,ii},1),ns);
   if jj == stimbefnum, c = cols(1,:); elseif jj==stimafternum, c=cols(2,:); else, c = 0.65*[1,1,1]; end
   plot((allspks{jj,ii}(id) + stimstartend(jj,1))/60, allamps{jj,ii}(id),'.','color',c);
   hold on;
   xline(stimstartend(jj,1)/60,'--r',jj);
end

yAx = [floor(min(cellfun(@min,allamps(:,ii)))/5)*5, ceil( max(cellfun(@max,allamps(:,ii)))/2)*2];
axis([0 expdur yAx]);
xticks(0:60:1000);               yticks(yAx(1):10:100);         box off;
%xlabel('Time (min)');   
ylabel('Std from noise');
title('All stimuli');
set(gca,'Xcolor','none');

subplot_tight(5,1,4,tt)
plot(expbins(1:end-1), binnedallspks,'color','k','linewidth',2);      
yAx = ceil(max(binnedallspks)/10)*10;
for jj = 1:size(allamps,1), xline(stimstartend(jj,1)/60,'--r',jj); end
axis([0 expdur 0 yAx]);
xticks(sort([0:60:1000,expdur]));               yticks(0:500:1e5);         box off;
xlabel('Time (min)');       
ylabel('spikes per min'); 


subplot_tight(5,3,[13 14],tt)

end



end