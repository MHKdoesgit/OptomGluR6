

function expinfo = loadExpInformation(datapath, varargin)
%
%%% loadExpInformation %%%
%
%
% This function opens a gui for getting information about the experiment. 
% It returns expinfo structure which contains all the necessary information
% about the details of the experiment. For more on the gui itself check the
% gui code in the appdesigner.
%
% ===============================Inputs====================================
%
%    datapath : path to experiment folder
%    expguiupdate : flag for updating an alreading existing file with the gui.
%
%================================Output====================================
%
%   expinfo : A structure containing parameters about the experiment time, 
%             setup, array, animal that was use and more.  
%
% written by Mohammad on 07.01.2021.

if nargin < 1
    datapath = uigetdir();
end

if nargin > 1,    expguiupdate = varargin{1}; else, expguiupdate = 0; end

expinfopath = [datapath,filesep, 'Data Analysis',filesep,'Raw Data'];
expinfofile = dir([expinfopath, filesep, 'experiment inforamtion for data recorded on *.mat']);
if ~isempty(expinfofile)
    expinfofile = [expinfopath, filesep, expinfofile.name];
end

if exist(expinfofile,'file') && ~expguiupdate 
    expinfo = load(expinfofile);
    return;
elseif exist(expinfofile,'file') && expguiupdate 
    expinfo = load(expinfofile);
    expinfo = updateexpgui(expinfo);
else
   expgui = experimentInfoGui;
   if ~isempty(datapath)
       expgui.Experimentpath.Value = datapath;
       expgui.UIFigure.UserData.datapath = datapath;
   end
   uiwait(expgui.UIFigure);
   expinfo = expgui.SaveButton.UserData;
   delete(expgui.UIFigure);
end

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function expinfo = updateexpgui(expinfo)
% open gui
expgui = experimentInfoGui;
% refill all the parts
expgui.Experimentpath.Value = expinfo.datapath;
expgui.UIFigure.UserData.datapath = expinfo.datapath;
expgui.kspath.Value = expinfo.kilosortpath;
expgui.frametimespath.Value = [fileparts(expinfo.kilosortpath),filesep,'frametimes'];
expgui.ExperimentDateDatePicker.Value = datetime (expinfo.expdate);
expgui.Project.Value = expinfo.project;
expgui.ExperimenterEF.Value = expinfo.experimenter;
expgui.CommentTextArea.Value = expinfo.comments;
radbut = expgui.StimulusParametersButtonGroup.Buttons;
expgui.StimulusParametersButtonGroup.SelectedObject = radbut (ismember({radbut.Text},expinfo.stimparatype));
expgui.StimulusParametersButtonGroup.SelectionChangedFcn(expgui, []);
expgui.SetupDropDown.Value = expinfo.setup.name;
expgui.SamplingRate.Value = expinfo.setup.samplingrate;
expgui.AmplifierDropDown.Value = expinfo.setup.amplifier;
expgui.Temperature.Value = expinfo.setup.temperature;
expgui.TemperatureGauge.Value = expinfo.setup.temperature;
expgui.FlowRateEditField.Value = expinfo.setup.flowrate;
% animal
expgui.AnimalKnob.Value = expinfo.animal.species;
expgui.AnimalKnob.ValueChangedFcn(expgui, []);
expgui.DateofBirthDatePicker.Value = datetime (expinfo.animal.dateofbirth);
if strcmpi(expinfo.animal.species,'marmoset')
    expgui.Genotype.Value = expinfo.animal.name;
    expgui.ID.Value = expinfo.animal.number;
else
    expgui.Genotype.Value = expinfo.animal.genotype;
    expgui.ID.Value = expinfo.animal.id;
end
expgui.Sex.Value = expinfo.animal.sex;
expgui.Age.Value = expinfo.animal.age;
% screen
expgui.ProjectorSwitch.Value = expinfo.screen.type;
expgui.RefreshRateKnob.Value = num2str (expinfo.screen.refreshrate);
if expinfo.screen.color, c='Color'; else, c= 'Black & White'; end
expgui.StimulatorSwitch.Value = c;
expgui.ScreenWidth.Value = expinfo.screen.resolution(1);
expgui.ScreenHeight.Value = expinfo.screen.resolution(2);
expgui.PixelSizeEF.Value = expinfo.screen.pixelsize;
expgui.BrightnessEditField.Value = num2str(expinfo.screen.brightness);
expgui.CalibrationDateDatePicker.Value = datetime (expinfo.screen.calibrationdate);
expgui.DelayEF.Value = expinfo.screen.delay;
% array
expgui.ArrayTypeDropDown.Value = expinfo.array.type;
expgui.NumberofElectrodes.Value = expinfo.array.numelectrodes;
expgui.ElectrodeDistance.Value = expinfo.array.electrodedistance;
expgui.ElectrodeSize.Value = expinfo.array.electrodesize;
expgui.ArrayName.Value = expinfo.array.name;
% eye
expgui.EyeSwitch.Value = expinfo.retina.eyenumber;
expgui.LeftRightEyeSwitch.Value = expinfo.retina.eyeside;
expgui.RecordingRegionKnob.Value = expinfo.retina.recordingregion;
% scores
expgui.BestScore.Value = expinfo.spikesorting.bestscore;
expgui.WorstScore.Value = expinfo.spikesorting.worstscore;
% options
if expinfo.options.parameterstxtfile, an = 'Yes'; else, an = 'No'; end
expgui.ParametersTextFileSwitch.Value = an;
if expinfo.options.excelfile, an = 'Yes'; else, an = 'No'; end
expgui.ExcelFileSwitch.Value = an;
expgui.SpikesorterSwitch.Value = expinfo.spikesorting.sorter;
expgui.FileTypeSwitch.Value = expinfo.spikesorting.datatype;
% stimulus names
stimnames = expinfo.stimulusnames;
ml = max(cellfun(@length,stimnames));
stimnameswarped = cellfun(@(x)strcat(x,repmat('.',1,ml-length(x))),stimnames,'un',0);
stimnameswarped = cellfun(@(x)([x(1:43),'...']),stimnameswarped,'un',0);
stimnameswarped = strrep(stimnameswarped,'.','');
expgui.StimulusListTextArea.Value = stimnameswarped;
expgui.UIFigure.UserData.stimulusnames = stimnames;
% wait for saving
uiwait(expgui.UIFigure);
expinfo = expgui.SaveButton.UserData;
delete(expgui.UIFigure);

end

