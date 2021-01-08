

function [ stimpara ] = readStimulusParameters( stimPath, screenWidth, screenHeight )
%READSTIMULUSPARAMETERS Creates a parameters struct based on a text file.
%==========================================================================
%   Input:
%       stimPath : path of the stimulus txt.
%       screenWidth : width of screen used during experiment.
%       screenHeight : Height of screen used during experiment.
%
%   Output:
%       stimpara: structure containing all stimulus parameters
%       expnum: experiment file data number or expId
%==========================================================================
%   Based on a function by Fernando and Mohammad.
%   Last updated 31.01.2017

[~,paramTxt]=fileparts(stimPath);
underInds = strfind(paramTxt, '_');
stimNumStr=paramTxt(1:underInds(1)-1);

fid = fopen(stimPath,'r');
C = textscan(fid, '%s','Delimiter','');
fclose(fid);
C = C{:};

stimpara.originalname = paramTxt(underInds(1)+1:end);
stimpara.expnumber = str2double(stimNumStr);

if numel(C)>0; desc=C{1}; else, desc='spontaneousactivity'; end
stimpara.stimulus = desc;
%get default parameters
stimpara = defaultStimulusParams( desc, stimpara, screenWidth, screenHeight );
%==========================================================================

%Overwrite all the defaults specified above, by reading the text file
idxEq=strfind(C, '=');
for ii=2:numel(C) 
    ide=idxEq{ii};
    if numel(ide)>0
        lineStr=C{ii};
        
        valueStr=strtrim(lineStr(ide+1:end));
        fieldStr=strtrim(lineStr(1:ide-1));

        if valueStr(end)==';'; valueStr=valueStr(1:end-1); end
          
        %take care of dots - FIX LATER!
        dp=strfind(fieldStr, '.');
        if numel( dp)==1
            ffieldstr=fieldStr(1:dp-1);
            sfieldstr=fieldStr(dp+1:end);
            
            %take care of booleans
            if strcmp(valueStr, 'true')
               stimpara.(ffieldstr).(sfieldstr)= true;
            elseif strcmp(valueStr, 'false')
                stimpara.(ffieldstr).(sfieldstr)= false;
            elseif numel(valueStr)>35 %take care of strings (e.g. paths)
                stimpara.(ffieldstr).(sfieldstr)= valueStr;
            else %then numeric value
                stimpara.(ffieldstr).(sfieldstr)=str2num(valueStr); %#ok
            end
             
             continue
         end

        %take care of booleans
        if strcmp(valueStr, 'true')
           stimpara.(fieldStr)= true;
        elseif strcmp(valueStr, 'false')
            stimpara.(fieldStr)= false;
        elseif numel(valueStr)>35 %take care of strings (e.g. paths)
            stimpara.(fieldStr)= valueStr;
        else %then numeric value
            stimpara.(fieldStr)=str2num(valueStr); %#ok
        end
        
        
        
    end   
end

%Check this one - DONE! IMPORTANT FOR CHECKERFLICKER( KEEP)
if ( isfield(stimpara,'stimulus') && ...
        any( strcmp({'pinknoise', 'checkerflicker', 'locallysparsenoise',...
        'FrozenNoise','locallysparsesubunitflash','checkerflickerplusmovie'},...
        stimpara.stimulus) ) )
    stimpara = applyScreen( stimpara, screenWidth, screenHeight );
end


%Change stimulus name according to parameters (different analysis)
%=================================================================
if strcmp({'checkerflicker'}, stimpara.stimulus) && stimpara.Nx*stimpara.Ny==1
    stimpara.stimulus='fullfieldflicker';
    stimpara.contrast=0.3;
end

if strcmp({'FrozenNoise'}, stimpara.stimulus) 
    if stimpara.Nx*stimpara.Ny==1
        stimpara.stimulus='frozenfullfieldflicker';
        stimpara.contrast=0.3;
    else
        stimpara.stimulus='frozencheckerflicker';
    end
end

if strcmp({'locallysparsenoise'}, stimpara.stimulus) && stimpara.stimduration<10
    stimpara.stimulus='locallysparsenoiseflicker';
end

%Long onoffsteps for iprgc
if strcmp({'onoffsteps'}, stimpara.stimulus) && stimpara.Nframes>600
    stimpara.stimulus='onoffstepslong';
end
if strcmp('gratingflashes', stimpara.stimulus)
    stimpara.stripewidths = str2num(stimpara.stripewidths); %#ok
    if numel(stimpara.Nphases) == 1
        stimpara.Nphases  = stimpara.Nphases * ones(size(stimpara.stripewidths));
    end
    if numel(stimpara.Norientations) == 1
        stimpara.Norientations  = stimpara.Norientations * ones(size(stimpara.stripewidths));
    end
end

%Movie for continuous image presentation
if strcmp({'imagesequence'}, stimpara.stimulus)
    spath=stimpara.path;
    stimstr='imagesequence';
    if (stimpara.trialduration==(stimpara.flashstop-stimpara.flashstart))
        stimstr='moviesequence';
    end
    if contains(spath, 'blur') || contains(spath, 'scales')
        stimpara.scales=getBlurScales(spath,stimstr);
        stimstr=[stimstr 'blur'];
    end
    stimpara.stimulus=stimstr;
end

%=================================================================

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [ newStimpara ] = applyScreen( stimpara, width, height )

newStimpara = stimpara;
newStimpara.Nx = ceil(width/stimpara.stixelwidth);
newStimpara.Ny = ceil(height/stimpara.stixelheight);
end
