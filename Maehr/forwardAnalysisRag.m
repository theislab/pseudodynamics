%% Preliminary
clear
close all;
clc;

%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Initial conditions (if not provided)
load dataMaehr8d_b_wN

options.grid = x;
optionsSim.grid = x;

% artificially shift cell at s=0 to s=e-14
D.ind.b1.hist{3}(153) = 10^(-20);

% compute mean u0 distribution
for i0 = 1:3
    u0(i0,:) = ksdensity(D.ind.b1.hist{i0},x,'support',[0,1],'function','pdf');
    N01(i0) = length(D.ind.b1.hist{i0});
    u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));%*D.ind.size(i0);
end
for i0 = 1:3
    %u02(i0,:) = ksdensity(D.ind.b2.hist{i0},x(20:end),'support',[0,1],'function','pdf');
    %N02(i0) = length(D.ind.b2.hist{i0});
    %u02(i0,:) = u02(i0,:)/trapz(x(20:end),u02(i0,:));%*D.ind.size(i0);
    u02(i0,:) = ksdensity(D.ind.b2.hist{i0},x(22:end),'support',[0,1],'function','pdf');
    N02(i0) = length(D.ind.b2.hist{i0});
    u02(i0,:) = u02(i0,:)/trapz(x(22:end),u02(i0,:));%*D.ind.size(i0);
end
w1=mean(N01./(N01+N02));
w2=mean(N02./(N01+N02));
optionsSim.u0 = [w1*mean(u0,1),w2*mean(u02,1)]*D.pop.mean(1);

% integral evaluted with trapezoidal rule * 1/h for averaging:
% (1/2*(u_i+u_{i+1})*h)*1/h
optionsSim.u0 = 0.5*([optionsSim.u0(1:n_grid-1)+optionsSim.u0(2:n_grid),optionsSim.u0(n_grid+1:end-1)+optionsSim.u0(n_grid+2:end)]);

%% Data
load dataMaehr_Rag_wN
D.pop.t = [0,2,4];
%for it = 3
%options.x_combined = union(options.grid(1:end-1)',D.Rag.xsd);
%D.Rag.csd_a = augment_cdf(D.Rag.xsd,options.x_combined,D.Rag.csd);
%end

%for it = 3
%    for ig = 1:n_grid-1
%        e_i = zeros(1,n_grid-1);
%        e_i(ig) = 1;
%        options.Aug_matrix(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined,e_i');
%    end
%end
D.Rag.csdr{2} = D.Rag.csdr{2}';
k1 = 2;
for it = 1:2
    options.x_combined{k1} = union(options.grid(1:end-1)',D.Rag.xsdr{it});
    D.csd_a{k1} = augment_cdf(D.Rag.xsdr{it},options.x_combined{k1},D.Rag.csdr{it});
                                 
    k1 = k1+1;
end
                                                                  
for it = 2:3
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{it},e_i');
    end
end
%% Definition of the Paramter Estimation Problem

% parameters
model = 'branching_froward';

% load estimation results
%Par =load('resultsMaehr');
Par =load('resultsRag');

for j=1:4
    %parameters = Par.Est{2}.parameters;
    parameters = Par.Est{j}.parameters;

    %modelfun = @ simulate_pd_branching_forward;
    modelfun = @(opt) computeSolForwardRag(parameters.MS.par(:,1),D,opt);

    alpha = [10,30,100,300];
    options.alpha = alpha(j);
    options.n_exp = 1;

    % Log-likelihood function
    objectiveFunction = @(opt) lsPseudodynamicsFA(modelfun(opt),D,options);

    %% evaluation of forward simulation
    k=1;
    for grid = 200:-1:1
        optionsSim.grid0(1) = 1;
        optionsSim.grid0(2) = grid;%29;%
        optionsSim.grid0(3) = grid;
        J(j,k) = objectiveFunction(optionsSim);
        k=k+1;
    end

    figure; plot(options.grid(100:1:299),J,'-o')
    ylabel('sum of squared residuals')
    xlabel('dpt')
end

save('fA_Rag')