% script for the compilation of *_syms files with AMICI
clear; close all; clc

%% COMPILATION

[exdir,~,~]=fileparts(which('pd_fv_wrap.m'));

% compile the model
tic;
amiwrap('pd_branching_fv','pd_branching_fv_syms',exdir)
t_wrap = toc

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/pd_branching_fv_syms']))