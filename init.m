warning('OFF', 'images:initSize:adjustingMag');

clear all;
close all;

convert_14_to_8bit('01-tif-16bit', '02-tif-8bit');
% FIXME it would be better to do this actually recursively and after reorganization
crop_recursive('02-tif-8bit', '03-cropped');
reorganize('03-cropped', '04-reorganized');
% mip('04-reorganized', '05-mips');
tophat('04-reorganized', '05-tophats');
mip('05-tophats', '06-mips');

mip('04-reorganized', '06-mips-full');

register('05-tophats', '06-mips', '07-registered');


register_cycles('04-reorganized', '06-mips', '07-registered-cycles');
%no: register('04-reorganized', '06-mips-full', '07-registered-full');
%register('04-reorganized', '05-mips', '06-registered');
%tophat('06-registered', '07-tophats');
%sequence('07-tophats', '08-sequence');
sequence('07-registered', '08-sequence');
filtered('08-sequence', '09-filtered');
centroids('09-filtered', '10-centroids');
%cp_comparison('07-tophats', '10-centroids', '11-cp-output', '12-comparison');
cp_comparison('07-registered', '10-centroids', '11-cp-output', '12-comparison');



sequence('XX-registered-gradient', '08-sequence');
filtered('08-sequence', '09-filtered');
centroids('09-filtered', '10-centroids');

cp_comparison('XX-registered-gradient', '10-centroids', '11-cp-output', '12-comparison');



register('04-reorganized', '06-mips-full', '07-registered');
tophat('07-registered', '08-tophats-registered');
%mip('08-tophats-registered', '09-mips-tophats');

sequence('08-tophats-registered', '09-sequence');
filtered('09-sequence', '07-registered', '10-filtered');
centroids('10-filtered', '11-centroids');

cp_comparison('07-registered', '09-mips-tophats', '11-centroids', '13-cp-output', '14-comparison');





