

function [ stimpara, expnum ] = loadStimulusParameters( filename, screenWidth, screenHeight )
%
%%% loadStimulusParameters %%%
%
%
% This function read the parameter file from experiments and extract the
% stiumuls properties used in the experiment, it works based on the named
% defied here.
%
%================================Inputs====================================
%
%    filename : stimulus name extracted from "parameters.txt".
%    screenWidth : width of screen used during experiment.
%    screenHeight : Height of screen used during experiment.
%
%================================Output====================================
%
%   stimpara : structure containing all stimulus parameters.
%   expnum: experiment file data number or expId
%
% Copied from loadStimdesc function from Fernando (2011-11-29),
% modified by Mohammad, 30.07.2014.
% updated for new stimuli by Mohammad 18.05.2015.
% updated to new naming system by Mohammad 13.07.2015.
% added new stimuli set on 13.11.2015.
% added new stimuli set on 27.12.2016.
% changed all sub-functions to private functions in a private folder on 28.12.2016.

expnum = 0;
[desc,numParams] = sscanf(filename, '%d_%s');

if (numParams < 2)
    desc = filename;
else
    expnum = desc(1);
    desc = char(desc(2:end))';
end

stimpara.originalname = desc;
stimpara.expnumber = expnum;

switch lower(desc)
    %---------------------------------------Spontaneous Activity------------------------------------
    case {'spontaneous','spontaneous_gray','spontaneous_activity_gray','spontaneous_activity_grey'...
            'spontaneous_black','spontaneous_activity_black','spontanous_activity'...
            'opto_spontaneous_activity_white_bleaching'...
            'opto_spontaneous_activity_grey_bleaching'...
            'opto_spontaneous_activity_grey_ndf3'...
            'spontaneous_activity_grey_beforedrug'...
            'spontaneous_activity_grey_washin','spontaneous_activity_grey_washout'...
            'spontactivity_gray_long','spontactivity_gray_adaptation'...
            'spontaneous_light','spontactivity_light','spontact_light'...
            'spontact_gray','spontactivity_gray','spontactivitygray','after2_spontact_gray'...
            'spontaneousactivity_gray','goingtogray','steppingtogray','goingbacktogray'...
            'spontactivitygray_long','turning_oled_on','spontactlight_adapt'...
            'spontactivity_dark', 'spontactivity_dark1','dark','spontactdark_duringadapt'...
            'spontactivity_darkadapt','spontactivity_darkadapt1','spontactivity_darkadapt2'...
            'spontaneousactivity_dark','spontaneous_dark', 'turningmonitoron','spontaneousdark'...
            'spontactivity_long_bright43','spontactivity_bright55','spontactivity_bright67'...
            'spontactivity_bright79','spontactivity_bright91','spontaneousactivity_gray_washout'...
            'spontactivity_gray_washin','spontactivity_gray_washout','spontactivity_gray_after',...
            'spontactivity_gray_afterlongtime','spontaneousactivity_darkadapt','spontaneousactivity_darkadapt2'}
        stimpara = spontaneous(desc, stimpara);
        %---------------------------------------Stim Sequence---------------------------------------
    case 'stimsequence'
        stimpara.stimulus = 'stimsequence';
        
        %---------------------------------------ON-OFF Steps----------------------------------------
    case{'onoffsteps30'
            'onoffsteps_contrast1_30frames'
            'onoffsteps30_50'
            'onoffstep30_50'
            'onoffsteps_60blinks_60preframes'
            'onoffsteps_60blinks'
            'onoffsteps30_contrast0.3'
            'onoffsteps30_preframes90_contrast0.4'
            'onoffsteps10_preframes30_contrast0.5'
            'onoffsteps10_preframes30_0.5contrast'
            'onoffstep30_50contrast'
            'onoffsteps30_50contrast'
            'onoffsteps30_preframes90_contrast1'
            'onoffsteps30_preframes90_contrast1_lap4'
            'onoffsteps30_preframes90_contrast1_cnqx'
            'onoffsteps30_preframes90_contrast1_washout'
            'onoffsteps_contrast1_30frames120preframes'
            'onoffsteps30_preframes120_contrast0.6'
            'onoffsteps30_preframes120_contrast1'
            'onoffsteps_100contrast_30frames_90pfr'
            'onoffsteps_100contrast_60frames_90pfr'
            'onoffsteps_100contrast_12frames_288pfr'
            'opto_onoffsteps_100contrast_60frames_90pfr'
            'onoffsteps_3600_for_iprgcs'
            'onoffsteps3600_for_iprgc'
            'opto_onoffsteps3600_for_iprgc'
            'onoffsteps1800_for_iprgc'
            'onoffsteps30preframes90contrast1p'
            'onoffsteps30preframes90contrasts1'
            'onoffsteps30preframes90contrast1'
            'onoffsteps3w0preframes90contrast1p'
            'onoffsteps30preframes60contrast1p'
            'onoffsteps30_preframes90_contrast1p'
            'onoffsteps30preframes90contrast1p1'
            'onoffsteps30preframes90contrast1_before1h'
            'onoffstep30preframes90contrast1_bright79'
            'onoffsteps30preframes90contrast1_bright91'
            'onoffsteps30preframes90contrast1_control'
            'onoffsteps30preframes90contrast1_am251_'
            'onoffstep30preframes90contrast1_after'
            'onoffsteps30preframes90contrast1_win55212'
            'onoffsteps30preframes90contrast1_after'
            'onoffsteps1200_preframes600_contrast0p6_iprgcs'
            'onoffsteps75_preframe150_contrast1_for75hz'
            'onoffsteps75_preframes150_contrast1_for75hz'}
        stimpara = loadDefault('onoffsteps', stimpara);
        stimpara = stimONOFFsteps(desc, stimpara);
        
        %------------------------------Full Field Flicker (Gaussian)--------------------------------
    case {'fff2blinks'
            'fff2bl'
            'fff_gauss_2blinks'
            'fff_gauss_2blinks_2nd'
            'fff_gauss_2blinks_3rd'
            'checkerflicker1600x1200_ga_2blinks'
            'ff_checkerflicker1600x1200_ga_2blinks'
            'checkerflicker1600x1200gauss2blinks'
            'fullfiledflicker_1600x1200_gauss_4blinks'
            'fullfieldflicker_1600x1200_gauss_1blink'
            'opto_fullfieldflicker_1600x1200_gauss_1blink'
            'fff1blink'
            'fff1blinks'
            'fff_gauss1blink'
            'fffgauss1blink'
            'fff1blinkafter'
            'fff1blink_duringadapt'
            'fff_gauss1blink_bright79'
            'fff1blink_bright91'}
        stimpara = loadDefault('fullfieldflicker', stimpara);
        stimpara = stimFullFieldFlicker(desc, stimpara);
        
        %-----------------------------------Direction Grating---------------------------------------
    case {'directiongratingslow'
            'directiongrating_slow'
            'directiongratingfast'
            'directiongrating_fast'
            'direction_grating_bar40_t80_10x5'
            'directiongrating_bar40_t80_10x5'
            'directiongrating_cyc10_bar300_dur400by100_bluegreen'
            'directiongrating_cyc10_bar300_dur400by100_cyanblack'
            'directiongrating_cyc10_bar300_dur400by100_bluebalck'
            'directiongrating_cyc10_bar300_dur400by100_greenblack'
            'direction_grating_fast_2hz_for_60hz'
            'direction_grating_slow_1hz_for_60hz'
            'directiongrating30_4hz'
            'directiongrating60_4hz'
            'directiongratingsequence_bw_6stimuli'
            'directiongratingsequence_color_4stimuli'
            'directiongratingsequence_color_3stimuli'
            'directiongratingsequence_bw_3stimuli_1contrast'
            'directiongratingsequence_color_3stimuli_12angles'
            'directiongratingsequence_color_3stimuli_8812angles'
            'directiongratingsequence_bw_2_stimuli_88angles'
            'opto_directiongratingsequence_bw_2_stimuli_88angles'
            'direction_grating32_4hz'
            'direction_grating60_4hz'
            'direction_grating60_2hz'
            'movingbars8dirs_bkg0_speed900'
            'directiongratingsequence_marmoset_6stimuli'
            'directiongratingsequence_marmoset_6stim'
            'directiongratingseq_marmoset_6stimuli'
            'directiongratingsequence_bw_6_stimuli_marmoset'
            'dssequence_2to5hz'
            'directiongratingsequence_marmoset_6stimuli_2to5hz'
            'directiongratingsequence_marmoset_6stimuli_2to5hz1'
            'directiongratingsequence_marmoset_6stimuli_2to5hz_after'
            'directiongratingseq_marmoset_6stimuli_2to5hz'
            'directiongrating_marmoset_6stimuli_2to5hz'
            'plaidpatches8dir_150deg_period80_pxl_v1_'
            'plaidpatches8dir_150deg_period80_pxl_v1'
            'directiongratingsequence_2stimuli_ds_os_for75hz'}
        stimpara = loadDefault('directiongrating', stimpara);
        stimpara = stimDirectionGrating(desc, stimpara);
        
        %-----------------------------White Noise (black and white)---------------------------------
    case {'checkerflicker10x10bw2blinks'
            'checkerflicker10x10bw2blinks2'
            'checkerflicker10x10bw_2blinks'
            'checkerflicker10x0bw_2blinks'
            'checkerflickerbw10x10_2bl'
            'checkerflicker10x10_bw2blinks'
            'checkerflicker8x8_bw2blinks'
            'checkerflicker5x5bw2blinks'
            'checkerflicker5x5_2blinks'
            'checkerflicker5x51blink'
            'checkerflicker8x8bw2blinks'
            'checkerflicker8x8bw2blincks'
            'checkerflicker1x1bw2blinks'
            'checkerflicker_8x8_bw_2blinks'
            'checkerflicker_10x10_bw_2blinks'
            'checkerflicker_8x8_bw_4blinks'
            'checkerflicker_3x3_bw_6blinks'
            'checkerflicker_2x2bw_4blinks'
            'checkerflicker_3x3bw_4blinks'
            'checkerflicker_8x8_bw_1blinks'
            'checkerflicker_8x8_bw_1contrast_1blink'
            'opto_checkerflicker_8x8_bw_1contrast_1blink'
            'checkerflickergabazine_10x10_bw_2blinks'
            'checkerflickertpmpa_10x10_bw_2blinks'
            'checkerflickerlap4_10x10_bw_2blinks'
            'checkerflickerstrychnine_10x10_bw_2blinks'
            'checkerflickerwashout_10x10_bw_2blinks'}
        stimpara = loadDefault('checkerflicker', stimpara);
        stimpara = stimCheckerflicker(desc, stimpara);
        
        %------------------------------stimCheckerFlickerPlusMovie----------------------------------
    case {'checkerflickerplusmovie2x2bw2blinks_doves60hz'
            'checkerflickerplusmovie2x2bw3blinks_doves60hz'
            'checkerflickerplusmovie4x4bw2blinks_doves60hz'
            'checkerflicker2x2bw2blinks_doves60hz1'
            'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed1'
            'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed2'
            'checkerflickerplusmovie2x2bw3blinks_doves60hz_seed3'}
        stimpara = stimCheckerFlickerPlusMovie(desc, stimpara);
        
        %------------------------------------Strips Flicker-----------------------------------------
    case{'whitenoisestripes'
            'whitenoisestripes_5x600'
            'stripeflicker_4x600_bw_4blinks'
            'stripeflicker2x600bw1blink_control'
            'stripeflicker2x600bw1blink_am251_'
            'stripeflicker2x600bw1blink_after'
            'stripeflicker2x600bw1blink_win55212'
            'stripeflicker2x600bw1blink'
            'stripeflicker5x600bw1blink'}
        stimpara = loadDefault('stripsflicker', stimpara);
        stimpara = stimStripsFlicker(desc, stimpara);
        
        %--------------------------.-Pink Noise & Temporal Pink Noise-------------------------------
    case {  'pinknoisestripes_0.1white'
            'pinkstripes0.1white'
            'pinknoisestripes_0.5white'
            'pinkstripes0.5white'
            'pinknoise10x10_2blinks_0.1white'
            'pinknoise10x10_2bl'
            'pinknoise10x10_2blinks'
            'pinknoise10x10_0.1white'
            'pinknoise0.1white'
            'temporalpinknoise'}
        stimpara = loadDefault('pinknoise', stimpara);
        stimpara = stimPinkStripes(desc, stimpara);
        
        %-------------------------------------Sub-unit Flashes--------------------------------------
    case {'subunitflash588_7lin'
            'subunitflash7_588_7lin'
            'subunitflash7_588lin24x10p0s2'
            'subunitflash21_196lin24x10p0s2'
            'subunitflash7_588log5p1s2'
            'subunitflash14_588lin20x12p0s2'
            'subunitflash7_588lin20x12p0s2'
            'suf_5pxspace_0.03_0.63_6cont_15x15'
            'suf_5pxspace_0.03_0.63_6cont_16angles_15x15'
            'suf_5pxspace_0.03_0.63_12cont_15x15'
            'suf_2pxspace_0.03_0.63_12con_16ang_15x15'
            'suf_nospace_0.03_0.63_12contrasts_16angles_15x15'
            'suf_5pxspace_0.03_0.63_6cont_7x7'
            'suf_2pxspace_0.03_0.63_12con_16ang_7x7'
            'suf_nospace_0.03_0.63_6cont_7x7'
            'suf_nospace_0.03_0.63_12contrasts_7x7'
            'suf_10px_0.03_0.63_12con_16ang_10x10'
            'suf_8px_0.03_0.63_12con_16ang_8x8'
            'suf_8px_0.03_0.63_12con_16ang_5x5'
            'subunitflash14x14ang24con10rep4dur200ms'
            'subunitflash14x14ang24con0rep4dur200ms'
            'subunitflash14x14ang24con10rep4dur200ms_500ms'}
        stimpara = stimSubUnitFlash(desc, stimpara);
        
        %--------------------------------------Image Sequences---------------------------------------
    case {'imgseq_512x512_300im_10r_mixed'
            'imgseq512x512_300im_10r_mixed'
            'imgseq512x512_40im_10r_scales5_15to240um'
            'imgseq512x512_40im_10r_scale5_15to240p'
            'imgseq512x512_40im_10r_scale5_15_to_240'
            'imgseq512x512_40im_10r_blur_scales5_30to180'
            'imgseq_movie800x600_doves30hz_blur_scales5'
            'imgseq_movies800x600_doves30hz_blur_scales5_15_to240_20r'}
        stimpara = stimImageSequence(desc, stimpara);
        
        %---------------------------------------Linear Ramps----------------------------------------
    case{'linearramps_200period'
            'linearramps_period500'
            'linearramp_period600'
            'linearramps_period900'
            'linearramps500period'
            'linearramps_period500frames'
            'linearramps600periodafterevenlonger'}
        stimpara.stimulus = 'linearramps';
        stimpara = stimLinearRamp (desc, stimpara);
        
        %------------------------------------Drifting Objects---------------------------------------
    case 'driftingobjects'
        stimpara.stimulus = 'linearramps';
        stimpara = stimDriftingObjects (desc, stimpara);
        
        %-------------------------------------Saccade Gratings--------------------------------------
    case {'saccadegrating_18pixshift2_100ms_for60hz'
            'saccadegrating_18pix_400fixshift2_100ms_for60hz'
            'saccadegratings_18pixshift2_100ms'
            'saccadegrating_16pixshift2_100ms'
            'saccadegrating_12pixshift2_67ms'
            'opto_saccadegrating_16pixshift2_100ms'
            'saccadegrating12pixshift2_50ms'
            'saccadegrating_12pixshift2_50'
            'saccadegrating_12pixshift2_50ms'
            'saccadegratings_12pixshift2_50ms'}
        stimpara = loadDefault('saccadegrating', stimpara);
        stimpara = stimSaccadeGrating (desc, stimpara);
        
        %-------------------------------------ON-OFF Grating----------------------------------------
    case {'onoffgrating60_for60hz'
            'onoffgrating60'
            'onoffgrating60d'
            'onoff30grating60d'
            'onoffgratingp60s30ph30ang0bluegreen'
            'onoffgratingp60s30ph30ang0cyanblack'
            'onoffgratingp60s30ph30ang0blueblack'
            'onoffgratingp60s30ph30ang0greenblack'
            'onoffgratingp32s30ph30ang0cyanblack'}
        stimpara = loadDefault('onoffgrating', stimpara);
        stimpara = stimOnOffGrating (desc, stimpara);
        
        %----------------------------Contrat Adaptation Full Field----------------------------------
    case {'ff_ca_ga_32_08_1blinks_sf2400'
            'ff_ca_ga_32_08_1blink_2400frames'
            'ff_ca_ga_32_08_2blinks_sf1200'
            'ff_ca_ga_32_08_2blinks'
            'ff_ca_ga_32_08_2blink_2400frames'
            'ff_ca_ga_32_04_1blinks_sf2400'
            'ff_ca_ga_32_04_1blink_2400frames'
            'ff_ca_ga_32_04_2blinks_sf1200'
            'ff_ca_ga_32_04_2blink_2400frames'
            'ff_ca_ga_32_16_1blinks_sf2400'
            'ff_ca_ga_32_16_1blink_2400frames'
            'ff_ca_ga_32_16_2blinks_sf1200'
            'ff_ca_ga_32_16_1blinks_2400frames'
            'ff_ca_ga_32_04_1blinks_2400frames'
            'ff_ca_ga_32_16_2blinks'
            'ff_ca_ga_32_16_2blink_2400frames'
            'ff_ca_ga_32_04_2blinks'
            'ff_ca_ga_32_08_1blinks2400frame'
            'ff_ca_ga_32_08_1blinks2400frames'
            'ff_ca_ga_32_16_1blinks2400frames'
            'ff_ca_ga_32_04_1blinks2400frames'
            'ff_ca_32_12_2blinks_sf1200_2blinks'
            'ff_ca_ga_32_08_1blinks_2400sf'}
        stimpara = loadDefault('contrastadaptation', stimpara);
        stimpara = stimContrastAdaptation(desc, stimpara);
        
        %------------------------------Spatial Contrast Adaptation----------------------------------
    case{'ch_ca_bw_100_20_1blinks_sf2400'
            'ch_ca_bw_100_20_1blinks_2400frames'
            'ch_ca_8x8_bw_100_20_1blinks_sf2400'
            'ch_ca_bw_100_20_2blinks_sf1200'
            'ch_ca_bw_100_20_2blinks'
            'ch_ca_8x8_bw_100_20_2blinks_sf1200'
            'ch_ca_2x2_bw_100_20_1blinks_sf2400'
            'ch_ca_4x4_bw_100_20_1blinks_sf2400'
            'ch_ca_6x6_bw_100_20_1blinks_sf2400'
            'ch_ca_12x12_bw_100_20_1blinks_sf2400'
            'ch_ca_16x16_bw_100_20_1blinks_sf2400'
            'ch_ca_20x20_bw_100_20_1blinks_sf2400'
            'ch_ca_24x24_bw_100_20_1blinks_sf2400'
            'ch_ca_32x32_bw_100_20_1blinks_sf2400'
            'ch_ca_bw_100_20_1blink_2400frames'
            'ch_ca_4x4_bw_100_20_1blinks_2400frames'
            'ch_ca_16x16_bw_100_20_1blinks_2400frames'
            'ch_ca_8x8_bw_100_20_1blinks_2400frames'
            'ch_4x4_ca_bw_100_20_1blink_2400frames'
            'ch_16x16_ca_bw_100_20_1blink_2400frames'
            'ch_gabazine_ca_bw_100_20_1blink_2400frames'
            'ch_ca_6x6_bw_100_20_1blinks_2400frames'
            'ch_ca_12x12_bw_100_20_1blinks_2400frames'
            'ch_ca_2x2_bw_100_20_1blinks_2400frames'
            'ch_ca_32x32_bw_100_20_1blinks_2400frames'
            'ch_ca_24x24_bw_100_20_1blinks_2400frames'
            'ch_ca_20x20_bw_100_20_1blinks_2400frames'
            'tpmpa_ch_ca_bw_100_20_1blinks_2400frames'
            'after_tpmpa_ch_ca_bw_100_20_1blinks_2400frames'
            'ch_ca_bw_100_20_1blinks_2400sf'
            'ca_5x5bw_100_20_1blink_sf1200_control'
            'ca_5x5bw1blink_sf1200_am251_'
            'ca_5x5bw_100_20_1blink_sf1200'
            'ca_5x5bw_100_20_1blink_sf1200_win55212'
            'ca_5x5bw_100_20_1blink_sf1200_after'
            % Alternate contrast
            'chalternate_ca_bw_100_20_1blinks_sf2400'
            'chalternate_ca_bw_100_20_2blinks'
            'chalternate_ca_bw_100_20_1blink_2400frames'
            'chalternate_ca_8x8_bw_100_20_1blinks_sf2400'
            'chalternate_ca_bw_100_20_2blinks_sf1200'
            'chalternate_ca_8x8_bw_100_20_2blinks_sf1200'
            'chalternate_ca_2x2_bw_100_20_1blinks_sf2400'
            'chalternate_ca_4x4_bw_100_20_1blinks_sf2400'
            'chalternate_ca_6x6_bw_100_20_1blinks_sf2400'
            'chalternate_ca_12x12_bw_100_20_1blinks_sf2400'
            'chalternate_ca_16x16_bw_100_20_1blinks_sf2400'
            'chalternate_ca_20x20_bw_100_20_1blinks_sf2400'
            'chalternate_ca_24x24_bw_100_20_1blinks_sf2400'
            'chalternate_ca_32x32_bw_100_20_1blinks_sf2400'
            'chalternate_ca_bw_100_20_2blink'
            'chalternate_ca_4x4_bw_100_20_1blinks_2400frames'
            'chalternate_ca_16x16_bw_100_20_1blinks_2400frames'
            'chalternate_ca_4x4_bw_100_20_1blink_2400frames'
            'chalternate_ca_16x16_bw_100_20_1blink_2400frames'
            'chalternate_ca_8x8_bw_100_20_1blink_2400frames'
            'chalternate_ca_bw_100_20_1blinks_2400sf'}
        stimpara = loadDefault('spatialcontrastadaptation', stimpara);
        stimpara = stimSpatialContrastAdaptation(desc, stimpara);
        
        %--------------------------------------RGB Steps--------------------------------------------
    case  {'rgbsteps50_contrast1redblack'
            'rgbsteps30contrast50redblack'
            'rgbsteps30_contrast50green'
            'rgbsteps30contrast50greenblack'
            'rgbsteps30contrast1greenblack'
            'rgbsteps50_contrast1blueblack'
            'rgbsteps30contrast50blueblack'
            'rgbsteps30contrast1blueblack'
            'rgbsteps30contrast1cyanblack'
            'rgbsteps30_contrast50purple'
            'rgbsteps30contrast50purpleblack'
            'rgbstepsopponent50_contrast1redblue'
            'rgbsteps30contrast50redblueoppnent'
            'rgbsteps30contrast1bluegreenoppnent'
            'rgbsteps30contrast1bluemeanisoresponse'
            'rgbsteps30contrast1greenmeanisoresponse'
            'rgbsteps30contrast1allrandomized'
            'rgbsteps_allrandomized_0.7contrast_30frames'
            'rgbsteps_allrandomized_1contrast_30frames'
            'rgbsteps_cyanblack_0.7contrast_30frames'
            'rgbsteps_cyanblack_1contrast_30frames'
            'rgbsteps_greenblack_0.7greencontrast_30frames'
            'rgbsteps_greenblack_1greencontrast_30frames'
            'rgbsteps_greenmean_0.7greencontrast_30frames'
            'rgbsteps_greenmean_1greencontrast_30frames'
            'rgbsteps_opponent_0.7contrast_30frames'
            'rgbsteps_opponent_1contrast_30frames'
            'rgbsteps_uvblack_0.7uvcontrast_30frames'
            'rgbsteps_uvblack_1uvcontrast_30frames'
            'rgbsteps_uvmean_0.7uvcontrast_30frames'
            'rgbsteps_uvmean_1uvcontrast_30frames'
            'rgbsteps_cyanblack_35contrast_30frames90pfr'
            'rgbsteps_cyanblack_35contrast_60frames90pfr'
            'rgbsteps_greenblack_80contrast_60frames90pfr'
            'rgbsteps_uvblack_80contrast_60frames90pfr'}
        stimpara = loadDefault('rgbsteps', stimpara);
        stimpara = stimRGBsteps(desc, stimpara);
        
        %----------------------------------White Noise Color----------------------------------------
    case {'checkerflicker8x8bw2blinks color'
            'checkerflicker8x8bw2blinkscolor'
            'checkerflicker10x10bw2blinkscolor'
            'checkerflicker_10x10bw_2blinks'
            'checkerflicker8x8bw2blinksredonly'
            'checkerflicker8x8bw2blinksblueonly'
            'checkerflicker8x8bw2blinksblackwhite'
            'checkerflicker8x8bw2blinkscyanonly'
            'checkerflicker8x8bw2blinksuvonly'
            'checkerflicker8x8bw2blinksgreenonly'
            'checkerflicker6x6bw2blinkscolor'
            'checkerflicker3x3bw2blinkscolor'
            'checkerflicker2x2bw2blinkscolor'
            'checkerflicker1x1bw2blinkscolor'
            'checkerflicker6x6bw2blinkscolor2seed'
            'checkerflicker6x6gauss2blinkscolor'
            'checkerflicker6x6bw2blinkscyanonly'
            'checkerflicker6x6bw2blinksblueonly'
            'checkerflicker6x6bw2blinksgreenonly'
            'checkerflicker6x6bw2blinksuvonly'
            'checkerflicker_4x4_bw_color_0.7contrast_2blinks'
            'checkerflicker_4x4_bw_cyanonly_0.7contrast_2blinks'
            'checkerflicker_4x4_bw_greenonly_0.7contrast_2blinks'
            'checkerflicker_4x4_bw_uvonly_0.7contrast_2blinks'
            'checkerflicker_4x4_gauss_color_0.7contrast_2blinks'
            'checkerflicker_6x6_bw_color_0.7contrast_2blinks'
            'checkerflicker_6x6_bw_cyanonly_0.7contrast_2blinks'
            'checkerflicker_6x6_bw_greenonly_0.7contrast_2blinks'
            'checkerflicker_6x6_bw_uvonly_0.7contrast_2blinks'
            'checkerflicker_6x6_gauss_color_0.7contrast_2blinks'
            'checkerflicker_8x8_bw_color_0.7contrast_2blinks'
            'checkerflicker_8x8_bw_cyanonly_0.7contrast_2blinks'
            'checkerflicker_8x8_bw_greenonly_0.7contrast_2blinks'
            'checkerflicker_8x8_bw_uvonly_0.7contrast_2blinks'
            'checkerflicker_8x8_gauss_color_0.3contrast_2blinks'
            'checkerflicker_6x6_bw_color_0.7contrast_6blinks'
            'checkerflicker_2x2_bw_color_0.7contrast_2blinks'
            'checkerflicker_2x2_bw_color_0.7contrast_4blinks'
            'checkerflicker_8x8_bw_color_0.7contrast_1blink'
            'checkerflicker_4x4_bw_color_0.7contrast_2blink'
            'checkerflicker_6x6_bw_color_0.7contrast_1blink'
            'checkerflicker_5x5_bw_3colors_1contrast_1blink'}
        stimpara = loadDefault('checkerflicker', stimpara);
        stimpara = stimCheckerflicker(desc, stimpara);
        
        %--------------------------Full Field Flicker Color (Gaussian)------------------------------
    case {'checkerflicker1600x1200gauss1blinks color'
            'checkerflicker1600x1200gauss1blinkcolor'
            'checkerflicker1600x1200gauss1blinkscolor'
            'checkerflicker1600x1200gauss2blinkscolor'
            'checkerflicker1600x1200gauss2blinksredonly'
            'checkerflicker1600x1200gauss2blinksblueonly'
            'checkerflicker1600x1200gauss1blinksredonly'
            'checkerflicker1600x1200gauss1blinksblueonly'
            'checkerflicker1600x1200gauss2blinkscyanonly'
            'checkerflicker1600x1200gauss2blinksgreenonly'
            'checkerflicker1600x1200gauss2blinksuvonly'
            'fullfieldflicker_1600x1200_gauss_color_2blinks'
            'fullfieldflicker_1600x1200_gauss_cyanonly_2blinks'
            'fullfieldflicker_1600x1200_gauss_greenonly_2blinks'
            'fullfieldflicker_1600x1200_gauss_uvonly_2blinks'
            'fullfieldflicker_1600x1200_gauss_color_1blinks'}
        stimpara = loadDefault('fullfieldflicker', stimpara);
        stimpara = stimFullFieldFlicker(desc, stimpara);
        
        %---------------------------------Color Intensity Steps-------------------------------------
    case{'colorintensitysteps120bluebkg30redsteps8skipval'
            'colorintensitysteps120redbkg30bluesteps8skipval'
            'colorintensitysteps60bluebkg30redsteps8skipval'
            'colorintensitysteps60redbkg30bluesteps8skipval'
            'colorintensitysteps120uvbkg60greensteps16skipval'
            'colorintensitysteps120greenbkg60uvsteps16skipval'
            'colorintensitysteps120greebbkg60uvsteps16skipval'
            'colorintensitysteps_120uvbkg_60greensteps_16skipval'
            'colorintensitysteps_120greenbkg_60uvsteps_16skipval'
            'colorintensitysteps_120uvbkg_60greensteps_10skipval'
            'colorintensitysteps_120greenbkg_60uvsteps_10skipval'}
        stimpara = loadDefault('colorintensitysteps', stimpara);
        stimpara = stimColorIntensitySteps(desc, stimpara);
        
        %----------------------------------Rotating Strips------------------------------------------
    case {'rotatingstrips6x6bw2blinkscolor2seed'
            'rotatingstrips6x6bw2blinkscyanonly'
            'rotatingstrips6x6bw2blinksgreenonly'
            'rotatingstrips6x6bw2blinksuvonly'
            'rotating_stripes_4x600_2blinks'
            'rotating_stripes_8x600_2blinks'
            'rotatingstripes_4x600_4blinks'
            'rotatingstripes_4x600_1blinks'
            'rotating_stripes_4x600_1blinks'
            'rotatingstripes4x600_1blink'
            'rotatingstripes_2x600_1blinks'
            'rotatingstripes2x600_1blinks'
            'rotatingstripes_2x600_1blink'
            'rotating_stripies_2x600_1blinks'
            'rotating_stripes_2x600_1blinks'
            'rotatingstripes5bw1blink5min'
            'rotatingstrips6x6gauss2blinkscolor2seed'
            'rotatingstrips6x6bw2blinks2seeduvonly'
            'rotatingstrips6x6bw2blinks2seedgreenonly'
            'rotatingstrips6x6bw2blinks2seedcyanonly'
            'rotatingstripes_2x600_bw_120preframes_18000dur'
            'rotatingstripes_4x600_bw_120preframes_18000dur'
            'rotatingstripes_8x600_bw_120preframes_18000sdur'
            'rotatingstripes_4x600_bw_color_120preframes_18000dur'
            'rotatingstripes_4x600_bw_color_120preframes_36000dur'
            'rotatingstripes_4x600_bw_cyanonly_120preframes_18000dur'
            'rotatingstripes_4x600_bw_greenonly_120preframes_18000dur'
            'rotatingstripes_4x600_bw_uvonly_120preframes_18000dur'
            'rotatingstripes_4x600_bw_color_120pfr_18000dur4blinks'
            'rotatingstripes_4x600_bw_color_120pfr_36000dur4blinks'
            'rotatingstripes_6x600_bw_color_120preframes_18000dur'
            'rotatingstripes_6x600_bw_cyanonly_120preframes_18000dur'
            'rotatingstripes_6x600_bw_greenonly_120preframes_18000dur'
            'rotatingstripes_6x600_bw_uvonly_120preframes_18000dur'
            'rotatingstripes_8x600_bw_color_120preframes_18000dur'
            'rotatingstripes_8x600_bw_cyanonly_120preframes_18000dur'
            'rotatingstripes_8x600_bw_greenonly_120preframes_18000dur'
            'rotatingstripes_8x600_bw_uvonly_120preframes_18000dur'}
        stimpara = loadDefault('rotatingstrips', stimpara);
        stimpara = stimRotatingStrips(desc, stimpara);
        
        %---------------------------------Contrast Step Ladder--------------------------------------
    case{'contraststeps'
            'contraststepladder'
            'contraststepladder_12fr_180preframes_white'
            'contraststepladder_12fr_180preframes_blue'
            'contraststepladder_12fr_180preframes_green'
            'contraststepladder_12fr_180preframes_cyan'
            'contraststepladderhighercontrasts'
            'contraststepladder_12fr_108preframes_cyan'
            'contraststepladder_12fr_108preframes_highcont'
            'contraststepladder12preframes48levels8'
            'contraststepladder_12fr_108pfr_1234681624325064contrasts'
            'contraststepladder_12fr_288pfr_2481624325064contrasts'
            'opto_contraststepladder_12fr_108pfr_1234681624325064c'
            'contraststepladder12preframe48levels8'
            'contraststepladder12preframes48levels8n'
            'cointraststepladder12preframes48levels8n'
            'contraststepsladder12preframes48levels8'
            'contraststepladder12preframes48_8levels'
            'contraststepladder12preframes48levels8_control'
            'contraststepsladder12preframes48levels8_am251'
            'contraststepladder12preframes48levels8_after'
            'contraststepladder12preframes48levels8_win55212'
            'contraststepladdere12preframes48levels8_after'
            'constraststepladder12preframes48levels8'
            'contrastladder_30blinks_60preframes'}
        stimpara = loadDefault('contraststepladder', stimpara);
        stimpara = stimContrastStepLadder(desc, stimpara);
        
        %-----------------------------------Silent Exchange-----------------------------------------
    case{'silentexchange_1_20g0uv_30frames120preframes'
            'silentexchange_2_18g-2uv_30frames120preframes'
            'silentexchange_3_16g-4uv_30frames120preframes'
            'silentexchange_4_14g-6uv_30frames120preframes'
            'silentexchange_5_12g-8uv_30frames120preframes'
            'silentexchange_6_20uv0g_30frames120preframes'
            'silentexchange_7_18uv-2g_30frames120preframes'
            'silentexchange_8_16uv-4g_30frames120preframes'
            'silentexchange_9_14uv-6g_30frames120preframes'
            'silentexchange_10_12uv-8g_30frames120preframes'}
        stimpara = loadDefault('silentexchange', stimpara);
        stimpara = stimSilentExchange(desc, stimpara);
        
        %---------------------------------Cone-Isolation Test---------------------------------------
    case{'coneisotest_1_0gmg0.5_30frames120preframes'
            'coneisotest_2_10gmg0.5_30frames120preframes'
            'coneisotest_3_20gmg0.5_30frames120preframes'
            'coneisotest_4_30gmg0.5_30frames120preframes'
            'coneisotest_5_40gmg0.5_30frames120preframes'
            'coneisotest_6_0gmg0.6_30frames120preframes'
            'coneisotest_7_10gmg0.6_30frames120preframes'
            'coneisotest_8_20gmg0.6_30frames120preframes'
            'coneisotest_9_30gmg0.6_30frames120preframes'
            'coneisotest_10_40gmg0.6_30frames120preframes'
            'coneisotest_10to40grcont_06newmean_70uvcont'}
        stimpara = loadDefault('coneisolationtest', stimpara);
        stimpara = stimConeIsolationTest(desc, stimpara);
        
        %---------------------------------Color Integration-----------------------------------------
    case{'colorintegration-20step2to20cont30stim120prefr'
            'colorintegration-20step2to20cont30stimdur120prefr'
            'colorintegration-20step2to20cont30stimdur90prefr'
            'colorintegration-40step4to40cont30stimdur90prefr'
            'sparsecolorintegration20x20gap2x2dur30circle'
            'sparsecolorintegration20x20gap2x2dur15circle'
            'sparsecolorintegration-40-4-40cont20x20gap2x2dur15circle'
            'colorintegrationbeforedrug-20-2-20cont30stim90prefr'
            'colorintegrationgabazine-20-2-20cont30stim90prefr'
            'colorintegrationtpmpa-20-2-20cont30stim90prefr'
            'colorintegrationlap4-20-2-20cont30stim90prefr'
            'colorintegrationstrychnine-20-2-20cont30stim90prefr'
            'colorintegrationwashout-20-2-20cont30stim90prefr'
            'colorintegration-20step2to20cont30stim90pfr_highphotopic'
            'colorintegration-60step6to60cont30stimdur90prefr'
            'colorintegration-60step2to60cont30stimdur90prefr'
            'colorintegration-20step2to20cont30stimdur90prefr_tpmpa'
            'colorintegrationhepes-20-2-20cont30stim90prefr'
            'colorintegrationhepes-20step2to20cont30stimdur90prefr'}
        stimpara = loadDefault('colorintegration', stimpara);
        stimpara = stimColorIntegration(desc, stimpara);
        
        %-----------------------------Color Integration Grating-------------------------------------
    case{'colorintegrationgrating60gw600stimdur60period'
            'colorintegrationgrating30gw600stimdur60period'
            'colorintegrationgrating6gw600stimdur60period'
            'colorintegrationgrating2gw600stimdur60period'
            'colorintegrationgrating60gw600-40-4-40stimdur60period'
            'colorintegrationgratinggabazine60gw600stimdur60period'
            'colorintegrationgratingtpmpa60gw600stimdur60period'
            'colorintegrationgratinglap460gw600stimdur60period'
            'colorintegrationgratingstrychnine60gw600stimdur60period'
            'colorintegrationgratingwashout60gw600stimdur60period'}
        stimpara = loadDefault('colorintegrationgrating', stimpara);
        stimpara = stimColorIntegrationGrating(desc, stimpara);
        
        %----------------------------------Color ISo-response---------------------------------------
    case{'chromisoresp36angle18cont3rep2to72cont30stim90pfr'
            'chromisoresp36angle10cont4rep2to42cont30stim90pfr'
            'chromisoresp36angle15cont4rep2to72cont15stim45pfr'
            'chromisoresp36angle15cont4rep2to72cont30stim90pfr'}
        stimpara = loadDefault('colorisoresponse', stimpara);
        stimpara = stimColorIsoresponse(desc, stimpara);
        
        %------------------------------Color ISo-response Grating-----------------------------------
    case{'chromisorg24ang9cont1rep2to72cont60gr600std60per'
            'chromisorg24ang7cont1rep2to42cont60gr600std60per'}
        stimpara = loadDefault('colorisoresponsegrating', stimpara);
        stimpara = stimColorIsoresponseGrating(desc, stimpara);
        
        %------------------------------------Color Exchange-----------------------------------------
    case{'colorexchange-40step20to40with20fix30stimdur90pfr'
            'colorexchange-40step10to40with10fix30stimdur90pfr'}
        stimpara = loadDefault('colorexchange', stimpara);
        stimpara = stimColorExchange(desc, stimpara);
        
        %------------------------------------Frozen Noise-------------------------------------------
    case {'frozennoise1600x1200gauss_rn1500_fn3002blinks'
            'frozennoise_1600x1200_gauss_4blinks_1200run900frozen'
            'checkerflicker_2x2bw4blinks_4200run_300frozen'
            'frozennoise1600x1200gauss_rn1500_fn300_4blinks'
            'frozennoise1600x1200gauss_rn1500_fn300_2blinks'
            'frozennoise1600x1200gauss_rn1500_fn300_1blink'
            'frozennoise1600x1200gauss_rgb_rn1500_fn300_1blink'
            'frozennoise1600x1200gauss_rn1500_fn3001blinks'
            'frozennoise1600x1200gauss_rn1500_fn300_1blinks'
            'frozennoise_1600x1200_gauss_rn1500_fn300_color_1blink'
            'fffgausswithrepeats_2blinks'
            'opto_frozennoise1600x1200gauss_rn1500_fn3001blinks'
            'opto_frozennoise1600x1200gauss_rn1500_fn3002blinks'
            'frozennoise8x8bw1blink1500run300freeze'
            'opto_frozennoise8x8bw1blink1500run300freeze'
            'frozennoise8x8bw2blinks1500run300freeze'
            'opto_frozennoise8x8bw2blinks1500run300freeze'
            'frozennoise5x5bw1blink1500run300freeze'
            'frozennoise5x5bwrgb1blink1500run300freeze'
            'frozennoise5x5bw1blinkrgb1500run300freeze'
            'frozennoisex5bw1blink1500run300freeze'
            'frozennoise5x5bw1blink_1500run_300freeze'
            'frozennoise5x5bw1blink1500runs300freeze'
            'frozennoise8x8binary_rn1500_fn300_1blinks'
            'frozennoise8x8bw1blink3000run600freeze'
            'frozennoise5x5bw1blink1500run300freeze_after'
            'checkerflicker9x9bw_3blinks_withrepeats'
            'frozennoise_8x8_bw_rn1500_fn300_color_1blink'
            'frozennoise5x5bw1blink3600run1200freeze'
            'frozennoise5x5bw1blink3600run1200freeze_40br'
            'frozennoise5x5bw1blink3600run1200freeze_0br'
            'fullfieldfrozennoise_gauss_1blink_1800run_600freeze'
            'fullfieldfrozennoise_gauss_1blink_1800run_600freeze_40'
            'fullfieldfrozennoise_gauss_1blink_1800run_600freeze_0'}
        stimpara = loadDefault('frozennoise', stimpara);
        stimpara = stimFrozenNoise(desc, stimpara);
        
        %-----------------------------------Chirp Stimulus------------------------------------------
    case {'chirpnf30pf120frq300lf05hf8ctr300frq5bluegreen'
            'chirpnf30pf120frq300lf05hf8ctr300frq5cyanblack'
            'chirpnf180pf240frq480lf0hf8ctr480frq2rep15euler'
            'opto_chirpnf180pf240frq480lf0hf8ctr480frq2rep15euler'
            'chirpstimulus'
            'chirpstimulus_lap4'
            'chirpstimulus_cnqx'
            'chirpstimulus_washout'
            'chirpstimulus_marmoset_15repeat'
            'chirpstimulus_marmoset_15'
            'chirpstimulus_marmoset_15r'
            'chirpstim_marmoset_15r'
            'chirpsimulus_marmoset_15r'
            'chirp_marmoset_15r'
            'chirpstimulus_15rp'
            'chirpstimulus_marmoset_15r_before1h'
            'chirpstimulus_marmoset_15r_bright79'
            'chirpstimulus_marmoset_15r_bright91'
            'chirpstimulus_marmoset_15r_after'
            'chirpstimulus_marmoset_15r_win55212'
            'chirpstimulus_bkg0contrast1_75hz'
            'chirpstimulus_bkg0contrast1_75hz_40'
            'chirpstimulus_bkg0contrast1_75hz_0'}
        stimpara = loadDefault('chirp', stimpara);
        stimpara = stimChirp(desc, stimpara);
        
        %----------------------------Temporal Frequency Sensitivity---------------------------------
    case {'tempfrequsensitivity_contrast0p8'
            'tempfrequsensitivity_contrast0p1'
            'tempfrequsensitivity_contrast0p8_lap4'
            'tempfrequsensitivity_contrast0p8_cnqx'
            'tempfrequsensitivity_contrast0p8_washout'
            'tempfrequsensitivity120preframes60levels8n'
            'tempfrequsensitivity120preframes60levels8'
            'tempfrequesensivity120preframes60levels8'
            'tempfrequencysensitivity120preframes60levels8n'
            'tempfrequsensitivity120preframes60level8n'
            'tempfreqsensitivity120preframes60levels8n'
            'tempfreqsensitivity120preframes60levels8'}
        stimpara = loadDefault('tempfreqsensitivity', stimpara);
        stimpara = stimTempFreqSensitivity(desc, stimpara);
        
        %--------------------Reversing Grating with Varying Spatial Period--------------------------
    case {'reversinggratingwithvaryingspatialperiod'
            'reversinggratingwvsp_30reversals_6widths'
            'reversinggratingwvsp12preframes60'
            'reversinggratingwvsp_30reversals_24816321000_6widths'
            'opto_reversinggratingwvsp_30rev_24816321000_6widths'
            'reversinggratingwvsp30preframes'
            'reversinggratingwvsp_12stim_25reversals_8widths'
            'reversinggratingwvsp12preframes60n'
            'reversinggratingwvp12preframes60'
            'reversingratingwvsp12preframes60'
            'reversingratingwvsp12preframes60f'
            'reversinggratingwvsp12preframes60f'
            'reversinggratingwvsp12preframes'
            'reversinggratingwvsp45preframes90_for75hz'
            'reversinggratingwvsp45preframes90_for75hz_40br'
            'reversinggratingwvsp45preframes90_for75hz_0br'}
        stimpara = stimSpatialReverseGrating(desc, stimpara);
        
        %--------------------------------Light Steps from Darkness----------------------------------
    case {'lightstepsfromdarkness60pfr240with10steps'
            'lightstepsfromdarkness60pfr240with10steps_washout'
            'lightstepsfromdarkness_60stim_240pfr_20steps_4repeats'
            'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats'
            'opto_lightstepsfromdarkness_60stim_240pfr_20steps_4rep'
            'lightstepsfromdarkness_75stim_300pfr_10steps_3repeats'
            'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats_40'
            'lightstepsfromdarkness_60stim_240pfr_20steps_3repeats_0'}
        stimpara = stimLightStepsFromDarkness(desc, stimpara);
        
        %----------------------------------Locally Sparse stimuli set ------------------------------
    case {'locallysparsenoise20x20gap2x2dur30circle'
            'locallysparsenoise20x20gap2x2dur60circle'
            'opto_locallysparsenoise20x20gap2x2dur60circle'
            'locallysparsenoise20x20gap2x2surround40dur2annulus'
            'locallysparsenoise20x20gap1x1surround40dur60annulus'
            'locallysparsenoise20x20gap2dur30'
            'locallysparsenoise20x20gap2dur30frames'
            'locallysparsenoise20x20gap2dur30f'
            'locallysparsenoise_20x20gap2dur30f'
            'locallysparsenoise_20x20gap2dur30'
            'locallysparsenoise20x20_gap2dur30'
            'locallysparsenoise20x20gap2dur30_am251_'
            'locallysparsenoise20x20gap2dur30f_control'
            'locallysparsenoise20x20_gap2_dur30_win55212'
            'locallysparsenoise20x20gap2dur30_after'
            'locallysparsenoise_20x20gap2dur30_after'
            'coloropponentcentersurround20cent60surr60stimdur'
            'coloropponentcentersurround20cent60surr60stimdur3x3gap'
            'coloropponentcentersurround25cent75surr30stimdur4x4gap'
            'coloropponentcentersurround20cent60surr15stimdur'
            'sparsebar50w6h4angle1gap60dur'
            'sparsespot8162550rad1gap60dur'
            'sparsespot8162550rad7gap60dur'
            'sparsespot162550rad1gap60dur'
            'sparsespot162550rad7gap60dur'
            'sparsespot162550rad4gap60dur'}
        stimpara = stimLocallySparseNoise(desc, stimpara);
        
        %---------SpatialFrequency & Contrast sensitivity with drifting grating stimuli set---------
    case {'spatialfrequencysensitivity'
            'contrastsensitivitywithdriftinggratings'
            'contrastsensitivitywithdriftinggrating'
            'contrastsensitivity_with_drifting_gratings'
            'contrastsensitivity_withdriftinggratings'
            'contrastsensitivity_with_driftinggratings'
            'contrastsensitivtywithdriftinggratings'
            'contrastsensitivtywithdriftinggrating'
            'contrastsensitivitywithdriftingratings_control'
            'contrastsensitivitywithdriftingratings_am251'
            'contrastsensitivitywithdriftinggratings_after'
            'contrastsensitivitywithdriftinggratings_control'
            'contrastsensitivitywithdriftinggratings_win55212'
            'contrastsensitivitydriftinggrating'
            'spatialfrequencysensitivity_bright55'
            'spatialfrequencysensitivity_bright67'
            'spatialfrequencysensitivity_bright79'
            'spatialfrequsensitivity_bright91'}
        stimpara = stimSpatialFrequencyContrastSensitivity(desc, stimpara);
        
        %----------------------------------fixation movie stimuli set-------------------------------
    case {'fixationmovie_agergb_part1', 'fixationmovie_agergb_part2','fixationmovie_agegray_part1',...
            'fixationmovie_agegray_part2', 'fixationmovie_mouse_for75hz375run1875freeze',...
            'fixationmovie_mouse_for75hz375run1875freeze_40br','fixationmovie_mouse_for75hz375run1875freeze_0br'}
        stimpara = stimFixationMovie(desc, stimpara);
        
        %----------------------------Omitted stimulus response stimuli set--------------------------
    case {'omittedstimulusresponse_black_on_white','omittedstimulusresponse_white_on_black'}
        stimpara = stimOmittedStimulusResponse(desc, stimpara);
        
        %----------------------------Objects moving background stimuli set--------------------------
    case {'omb_bg1x1corr4stdvc150_gsteps3stdv','omb_bg4x4corr4stdvc150_gsteps3stdv',...
            'omb_bg4x4corr8stdv_c150_gsteps3stdv',...
            'omb_bg1x1corr4stdv_c150_gsteps3stdv_15min','omb_bg4x4corr4stdv_c150_gsteps3stdv_15min'}
        stimpara = stimObjectsMovingBackground(desc, stimpara);
        
        %---------------------------------------Bar code stimuli set--------------------------------
    case {'barcodestimulus_plusvertical_120preframe'
            'barcodestimulus_seed20000_for75hz'}
        stimpara = stimBarCode(desc, stimpara);
        
    otherwise
        warning('StimDesc:FileNotFound',['Uh oh, I dont know that file: ' filename]);
end


if ( isfield(stimpara,'stimulus') && any( strcmp({'pinknoise', 'checkerflicker'}, stimpara.stimulus) ) )
    stimpara = applyScreen( stimpara, screenWidth, screenHeight );
end

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [ newStimpara ] = applyScreen( stimpara, width, height )

newStimpara = stimpara;
newStimpara.Nx = ceil(width/stimpara.stixelwidth);
newStimpara.Ny = ceil(height/stimpara.stixelheight);
end
