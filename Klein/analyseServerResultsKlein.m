%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Initial conditions (if not provided)
load dataKlein

optionsSim.grid = x;
options.grid = x;

% artificially shift cell at s=0 to s=e-14
D.ind.hist{1}(22) = 10^(-20);

% compute mean u0 distribution
options.u0 = ksdensity(D.ind.hist{1},x,'support',[0,1],'function','pdf');

options.u0 = 0.5*([options.u0(1:n_grid-1)+options.u0(2:n_grid),options.u0(n_grid+1:end-1)+options.u0(n_grid+2:end)]);

modelfun = @ simulate_pd_fv_Klein;

k1 = 1;
for it = 1:4
    options.x_combined{k1} = union(options.grid(1:end-1)',D.xsdt{it});
    D.csd_a{k1} = augment_cdf(D.xsdt{it},options.x_combined{k1},D.csdt{it});
                                                                
    k1 = k1+1;
end
                                                                
for it = 2:length(D.pop.t)
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{it},e_i');
    end
end

rep{1}=load('./parametersKlein/parameters_Klein');

for i = 1:12
[~,ind] = sort(rep{1}.parameters.MS.logPost((i-1)*40+1:i*40),'descend');
ind = ind([find(~isnan(rep{1}.parameters.MS.logPost(ind+40*(i-1)))),...
           find( isnan(rep{1}.parameters.MS.logPost(ind+40*(i-1))))]);

Est{i}.parameters.number = rep{1}.parameters.number;
Est{i}.parameters.min = rep{1}.parameters.min;
Est{i}.parameters.max = rep{1}.parameters.max;
Est{i}.parameters.constraints = rep{1}.parameters.constraints;
Est{i}.parameters.name = rep{1}.parameters.name;
Est{i}.parameters.guess = rep{1}.parameters.guess;

if isfield(rep{1}.parameters.MS,'par0')
    Est{i}.parameters.MS.par0 = rep{1}.parameters.MS.par0(:,40*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'par')
    Est{i}.parameters.MS.par = rep{1}.parameters.MS.par(:,40*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'logPost0')
    Est{i}.parameters.MS.logPost0 = rep{1}.parameters.MS.logPost0(40*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'logPost')
    Est{i}.parameters.MS.logPost = rep{1}.parameters.MS.logPost(40*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'gradient')
    Est{i}.parameters.MS.gradient = rep{1}.parameters.MS.gradient(:,40*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'hessian')
    Est{i}.parameters.MS.hessian = rep{1}.parameters.MS.hessian(:,:,40*(i-1)+ind);
end

if isfield(rep{1}.parameters.MS,'n_objfun')
    Est{i}.parameters.MS.n_objfun = rep{1}.parameters.MS.n_objfun(40*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'n_iter')
    Est{i}.parameters.MS.n_iter = rep{1}.parameters.MS.n_iter(40*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'t_cpu')
    Est{i}.parameters.MS.t_cpu = rep{1}.parameters.MS.t_cpu(40*(i-1)+ind)';
end

if isfield(rep{1}.parameters.MS,'exitflag')
    Est{i}.parameters.MS.exitflag = rep{1}.parameters.MS.exitflag(40*(i-1)+ind)';
end

sol{i,1} = modelfun(D.pop.t,Est{i}.parameters.MS.par(:,1),options.u0,[],[]);
sol{i,1}.p = bsxfun(@rdivide,sol{i,1}.y(:,1:299),sol{i,1}.y(:,end));

end
                                   
k=1
for i = [1,5,9]
    solGIF{k,1} = modelfun(0:0.1:7,Est{i}.parameters.MS.par(:,1),options.u0,[],[]);
    solGIF{k,1}.p = bsxfun(@rdivide,solGIF{k,1}.y(:,1:299),solGIF{k,1}.y(:,end));
    k=k+1;
end

for i=1:12
    alpha = [0,0,0,0,1,1,1,1,10,10,10,10];
    options.alpha = alpha(i);
    options.n_exp = 1;
    [negLogL(i)] = lsPseudodynamicsKlein(Est{i}.parameters.MS.par(:,1),modelfun,D,options);
end

A = reshape(negLogL,4,3);
