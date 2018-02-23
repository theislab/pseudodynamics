
n_grid1 = 300;
grid = linspace(0,1,n_grid1);

% (1/h)^2
h2inv = (1/(grid(2)-grid(1)))^2;

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

load('dataQuasiSteady10000')
load('resultsQuasiSteady')

%% parameters
for i=1:3
    par{i} = plotParameters(Est{i}.parameters.MS.par(:,1));
end

%% CI local approximation

load('resultsQuasiSteady_mS4')

optionsCI.mode = 'text';
alpha = [0.9,0.95,0.99];
for i=1:3
    Est{i}.parameters = getParameterConfidenceIntervals(Est{i}.parameters, alpha, optionsCI);
end
%%
clear

n_grid1 = 300;
grid = linspace(0,1,n_grid1);

% (1/h)^2
h2inv = (1/(grid(2)-grid(1)))^2;

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

load('dataAttractor10000')
load('resultsAttractor')
%% parameters
for i=1:3
    par{i} = plotParameters(Est{i}.parameters.MS.par(:,1));
end
%% CI local approximation
load('resultsAttractor')

alpha = [0.9,0.95,0.99];
optionsCI.mode = 'text';
for i=1:3
    Est{i}.parameters = getParameterConfidenceIntervals(Est{i}.parameters, alpha, optionsCI);
end