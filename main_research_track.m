clear all
close all
clc
%% constraints
WIN = 50; %window size for cross-correlogram (ms?)
DELTA = 1; % bin width for histogram (ms?)
NPAR = 102;
MAX = 200000;
%% plot histogram
T = 5400; %duration of recordin [s]
cc_list = linear_crossCorrelogram(cell4,cell9,T);
bin_width = DELTA;
bin_num = 2* WIN / bin_width;

histogram(cc_list{1},bin_num,"BinEdges",-w:w,'Normalization','count');
xlim([-w,w]);
xlabel('Spike Time Differences [ms]');
ylabel('Count');
title('Cross-Correlogram');

