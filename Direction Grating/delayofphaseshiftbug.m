

function stimstarts = delayofphaseshiftbug(para)
%--------------------------------------------------------------------------
%Get delays for correction
% Nframes = sum( para.numstimuli*para.preframestimuli+...
%     para.Nangles.*para.cycles.*(para.regeneration + para.duration));

relframecount = 0;
cyclecount = 0;
curmodule = 0;
curstimulus = 0;
k = 0;
allstarts = zeros(sum(para.Nangles.*para.cycles),1);
indstart = 1;
while true
    moduleDuration = para.duration(curstimulus+1) + para.regeneration(curstimulus+1);
    
    relmodulecount = relframecount - para.preframestimuli;
    modulePhase = mod(relmodulecount, moduleDuration);
    
    if modulePhase == para.regeneration(curstimulus+1)
        %drawGrating(dsvars, squareWave[curstimulus]);	// draw grating based on the parameters
        allstarts(indstart) = mod(k/para.period(curstimulus+1),1);
        indstart = indstart +1;
    end
    
    relframecount = relframecount + 1;
    if (modulePhase + 1) >= moduleDuration
        curmodule = curmodule + 1;
    end
    
    if curmodule >= para.Nangles(curstimulus+1)
        curmodule = 0;
        cyclecount = cyclecount + 1;
    end
    
    if cyclecount >= para.cycles(curstimulus+1)
        cyclecount = 0;
        relframecount = 0;
        curstimulus = curstimulus + 1;
    end
    k = k +1;
    if curstimulus >=para.numstimuli
        break;
    end
end

stimpulses = [0 cumsum(para.cycles .* para.Nangles)];

stimstarts = cell(para.numstimuli,1);
for istim = 1:para.numstimuli
    stimdelays = reshape(allstarts(stimpulses(istim)+1:stimpulses(istim+1)),...
        para.Nangles(istim),para.cycles(istim));
    stimstarts{istim} = stimdelays*para.period(istim)/para.fps;
end

end


