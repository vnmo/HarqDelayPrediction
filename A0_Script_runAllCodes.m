%% All codes for generating the figures can be run together with this script.
%% Comment/uncomment whichever is necesary
%% This could also be used as a starting point to find all functions


close all;
h1 = figure(1)
h2 = figure(2)
h3 = figure(3)

clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
bler_vs_dByRttRatio
% 
clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
blerVsNrb_vsDelayTarget
% 
clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
blerVsNrb_vsNrbPerPktLen

clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
blerVsNrb_vsPktLen

clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
blerVsNrb_vsRtt

% clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
% blerVsNrb_vsSnr

% clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
% blerVsPktLen_vsSnr

% clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
% blerVsRtt_vsSnr 

clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
throughputVsf_vsNrb
throughputvsf_vsRTT

clf(h1,'reset');clf(h2,'reset');clf(h3,'reset')
perf_comparison