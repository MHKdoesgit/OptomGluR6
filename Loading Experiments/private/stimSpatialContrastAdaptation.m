

function [ newStimpara ] = stimSpatialContrastAdaptation(desc, stimpara)
%
%%% stimSpatialContrastAdaptation %%%
%
%
% This function is a private function of loadStimulusParameters function.
%
%================================Inputs====================================
%
%   desc : original name of the mcd file and experiment data.
%   stimpara : stimulus parameter structure.
%
%================================Output====================================
%
%   newStimpara : updated stimpara structure.
%

switch lower(desc)
    
    %Spatial Contrast Adaptation
    case {'ch_ca_bw_100_20_1blinks_sf2400'
            'ch_ca_bw_100_20_1blinks_2400frames'
            'ch_ca_8x8_bw_100_20_1blinks_sf2400'
            'ch_ca_bw_100_20_1blink_2400frames'
            'ch_ca_8x8_bw_100_20_1blinks_2400frames'
            'ch_gabazine_ca_bw_100_20_1blink_2400frames'
            'tpmpa_ch_ca_bw_100_20_1blinks_2400frames'
            'after_tpmpa_ch_ca_bw_100_20_1blinks_2400frames'
            'ch_ca_bw_100_20_1blinks_2400sf'}
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_bw_100_20_2blinks_sf1200'
            'ch_ca_bw_100_20_2blinks'
            'ch_ca_8x8_bw_100_20_2blinks_sf1200'}
        stimpara.Nblinks = 2;
        stimpara.switchFrames = 1200;
        stimpara.twofixed = true;
        
    case {'ch_ca_2x2_bw_100_20_1blinks_sf2400'
            'ch_ca_2x2_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.graywidth = 2;
        stimpara.grayheight = 2;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_4x4_bw_100_20_1blinks_sf2400'
            'ch_ca_4x4_bw_100_20_1blinks_2400frames'
            'ch_4x4_ca_bw_100_20_1blink_2400frames'}
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
        stimpara.graywidth = 4;
        stimpara.grayheight = 4;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_6x6_bw_100_20_1blinks_sf2400'
            'ch_ca_6x6_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
        stimpara.graywidth = 6;
        stimpara.grayheight = 6;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_12x12_bw_100_20_1blinks_sf2400'
            'ch_ca_12x12_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 12;
        stimpara.stixelheight = 12;
        stimpara.graywidth = 12;
        stimpara.grayheight = 12;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_16x16_bw_100_20_1blinks_sf2400'
            'ch_ca_16x16_bw_100_20_1blinks_2400frames'
            'ch_16x16_ca_bw_100_20_1blink_2400frames'}
        stimpara.stixelwidth = 16;
        stimpara.stixelheight = 16;
        stimpara.graywidth = 16;
        stimpara.grayheight = 16;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_20x20_bw_100_20_1blinks_sf2400'
            'ch_ca_20x20_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 20;
        stimpara.stixelheight = 20;
        stimpara.graywidth = 20;
        stimpara.grayheight = 20;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_24x24_bw_100_20_1blinks_sf2400'
            'ch_ca_24x24_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 24;
        stimpara.stixelheight = 24;
        stimpara.graywidth = 24;
        stimpara.grayheight = 24;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
    case {'ch_ca_32x32_bw_100_20_1blinks_sf2400'
            'ch_ca_32x32_bw_100_20_1blinks_2400frames'}
        stimpara.stixelwidth = 32;
        stimpara.stixelheight = 32;
        stimpara.graywidth = 32;
        stimpara.grayheight = 32;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = true;
        
        
    case {'ca_5x5bw_100_20_1blink_sf1200_control'
            'ca_5x5bw1blink_sf1200_am251_'
            'ca_5x5bw_100_20_1blink_sf1200'
            'ca_5x5bw_100_20_1blink_sf1200_win55212'
            'ca_5x5bw_100_20_1blink_sf1200_after'}
        stimpara.stixelwidth = 5;
        stimpara.stixelheight = 5;
        stimpara.graywidth = 5;
        stimpara.grayheight = 5;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 1200;
        stimpara.twofixed = true;
        
        % Spatial Contrast Adaptation Alternate
    case {'chalternate_ca_bw_100_20_1blinks_sf2400'
            'chalternate_ca_bw_100_20_1blink_2400frames'
            'chalternate_ca_8x8_bw_100_20_1blinks_sf2400'
            'chalternate_ca_8x8_bw_100_20_1blink_2400frames'
            'chalternate_ca_bw_100_20_1blinks_2400sf'}
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case {'chalternate_ca_bw_100_20_2blinks_sf1200'
            'chalternate_ca_8x8_bw_100_20_2blinks_sf1200'
            'chalternate_ca_bw_100_20_2blink'
            'chalternate_ca_bw_100_20_2blinks'}
        stimpara.Nblinks = 2;
        stimpara.switchFrames = 1200;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_2x2_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 2;
        stimpara.stixelheight = 2;
        stimpara.graywidth = 2;
        stimpara.grayheight = 2;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case {'chalternate_ca_4x4_bw_100_20_1blinks_sf2400'
            'chalternate_ca_4x4_bw_100_20_1blinks_2400frames'
            'chalternate_ca_4x4_bw_100_20_1blink_2400frames'}
        stimpara.stixelwidth = 4;
        stimpara.stixelheight = 4;
        stimpara.graywidth = 4;
        stimpara.grayheight = 4;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_6x6_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 6;
        stimpara.stixelheight = 6;
        stimpara.graywidth = 6;
        stimpara.grayheight = 6;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_12x12_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 12;
        stimpara.stixelheight = 12;
        stimpara.graywidth = 12;
        stimpara.grayheight = 12;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case {'chalternate_ca_16x16_bw_100_20_1blinks_sf2400'
            'chalternate_ca_16x16_bw_100_20_1blinks_2400frames'
            'chalternate_ca_16x16_bw_100_20_1blink_2400frames'}
        stimpara.stixelwidth = 16;
        stimpara.stixelheight = 16;
        stimpara.graywidth = 16;
        stimpara.grayheight = 16;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_20x20_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 20;
        stimpara.stixelheight = 20;
        stimpara.graywidth = 20;
        stimpara.grayheight = 20;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_24x24_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 24;
        stimpara.stixelheight = 24;
        stimpara.graywidth = 24;
        stimpara.grayheight = 24;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    case 'chalternate_ca_32x32_bw_100_20_1blinks_sf2400'
        stimpara.stixelwidth = 32;
        stimpara.stixelheight = 32;
        stimpara.graywidth = 32;
        stimpara.grayheight = 32;
        stimpara.Nblinks = 1;
        stimpara.switchFrames = 2400;
        stimpara.twofixed = false;
        
    otherwise
        warning(['Uh oh, I dont know that file: ' desc])
end

newStimpara = stimpara;

end


