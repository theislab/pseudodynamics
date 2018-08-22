% analyze the results of the example
n_grid1 = 300;
grid = linspace(0,1,n_grid1);
% (1/h)^2
h2inv = (1/(grid(2)-grid(1)))^2;
grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

load('dataExample')
load('./results/resultsExample')
%% plot MultiStarts
for i=1:3
    plotMultiStarts(Est{i}.parameters)
end
%% plot parameters splines
for i=1:3
    par{i} = computeParameters(Est{i}.parameters.MS.par(:,1));
    
    figure; 
    plot(par{i}.grid,par{i}.D)
    xlabel('cell state')
    ylabel('diffusion parameter')
    legend('estimated')

    figure; 
    plot(par{i}.grid,par{i}.v)
    xlabel('cell state')
    ylabel('drift parameter')
    legend('estimated')

    figure; 
    plot(par{i}.gridx,par{i}.a)
    xlabel('cell state')
    ylabel('birth-death parameter')
    legend('estimated')
end

%% CI approximation using profile likelihoods
load('parameterProfilesExample')

alpha = [0.9,0.95,0.99];
parameters = getParameterConfidenceIntervals(parameters, alpha);
