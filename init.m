warning('OFF', 'images:initSize:adjustingMag');

clear all;
close all;

% OBS! this is not really intended to be run as a script
% run each step by itself and check the results

% convert 14-bit inputs to 8-bits
convert_14_to_8bit('01-tif-16bit', '02-tif-8bit');
% crop off the MIP-related artifacts
crop_recursive('02-tif-8bit', '03-cropped');
% reorganize for a consistent directory structure
reorganize('03-cropped', '04-reorganized');
% produce MIPs for the registration process
mip('04-reorganized', '05-mips');
% register the hybridization cycles
register('04-reorganized', '05-mips', '06-registered');
% produce top-hats of the registered base images
tophat('06-registered', '07-tophats-registered');
% perform sequencing
sequence('07-tophats-registered', '08-sequence');
% exclude pixels not matching the configured thresholds
filtered('08-sequence', '06-registered', '09-filtered');
% find the relevant centroids
centroids('09-filtered', '10-centroids');
%
% here you should run the CellProfiler pipeline and output to the directory 11-cp-output
%
% compare the results
cp_comparison('06-registered', '05-mips', '10-centroids', '11-cp-output', '12-comparison');
