% sort the collect individual starts
%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Initial conditions (if not provided)
load dataExample

optionsSim.grid = x;
options.grid = x;

% compute mean u0 distribution
for i0 = 1:3
    u0(i0,:) = ksdensity(D.ind.hist{i0},x,'support',[0,1],'function','pdf');
    N01(i0) = length(D.ind.hist{i0});
    u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));
end

options.u0 = mean(u0,1)*D.pop.mean(1);
options.u0 = 0.5*(options.u0(1:n_grid-1)+options.u0(2:n_grid));

options.n_exp = 1;

rep{1}=load('./results/parametersExample.mat');

%%
for i = 1:3
[~,ind] = sort(rep{1}.parameters.MS.logPost((i-1)*80+1:i*80),'descend');
ind = ind([find(~isnan(rep{1}.parameters.MS.logPost(ind+80*(i-1)))),...
           find( isnan(rep{1}.parameters.MS.logPost(ind+80*(i-1))))]);

Est{i}.parameters.number = rep{1}.parameters.number;
Est{i}.parameters.min = rep{1}.parameters.min;
Est{i}.parameters.max = rep{1}.parameters.max;
Est{i}.parameters.constraints = rep{1}.parameters.constraints;
Est{i}.parameters.name = rep{1}.parameters.name;
Est{i}.parameters.guess = rep{1}.parameters.guess;

if isfield(rep{1}.parameters.MS,'par0')
    Est{i}.parameters.MS.par0 = rep{1}.parameters.MS.par0(:,80*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'par')
    Est{i}.parameters.MS.par = rep{1}.parameters.MS.par(:,80*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'logPost0')
    Est{i}.parameters.MS.logPost0 = rep{1}.parameters.MS.logPost0(80*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'logPost')
    Est{i}.parameters.MS.logPost = rep{1}.parameters.MS.logPost(80*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'gradient')
    Est{i}.parameters.MS.gradient = rep{1}.parameters.MS.gradient(:,80*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'hessian')
    Est{i}.parameters.MS.hessian = rep{1}.parameters.MS.hessian(:,:,80*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'n_objfun')
    Est{i}.parameters.MS.n_objfun = rep{1}.parameters.MS.n_objfun(80*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'n_iter')
    Est{i}.parameters.MS.n_iter = rep{1}.parameters.MS.n_iter(80*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'t_cpu')
    Est{i}.parameters.MS.t_cpu = rep{1}.parameters.MS.t_cpu(80*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'exitflag')
    Est{i}.parameters.MS.exitflag = rep{1}.parameters.MS.exitflag(80*(i-1)+ind)';
end

    options_sim.sensi = 1;
    options_sim.maxsteps = 1e6;
    options_sim.sx0 = zeros(n_grid-1,length(rep{1}.parameters.MS.par(:,1)));
    sol{i} = simulate_pd_fv(D.pop.t,rep{1}.parameters.MS.par(:,(i-1)*80+ind(1)),options.u0,[],options_sim);
    sol{i}.p = bsxfun(@rdivide,sol{i}.y(:,1:299),sol{i}.y(:,end));
end

save('resultsExample','Est','sol')