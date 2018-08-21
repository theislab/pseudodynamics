% script for the compilation of *_syms files with AMICI
clear; close all; clc

%% COMPILATION

[exdir,~,~]=fileparts(which('compartment_wrap.m'));

% compile the model
tic;
amiwrap('compartment','compartment_syms',exdir)
t_wrap = toc

% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/compartment_syms']))