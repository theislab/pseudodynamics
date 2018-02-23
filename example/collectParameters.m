clear;
% put path to results of individual multi starts here
path = pwd;
par{1} = load([pwd '/parametersToy_' num2str(1) '.mat']);
parameters = par{1}.parameters;
options = par{1}.options;
optionsMultistart = par{1}.optionsMultistart;
k=2;
for i=2:240
    try
        par{i} = load([path '/parametersToy5_' num2str(i) '.mat']);
        parameters.MS.par(:,k) = par{i}.parameters.MS.par;
        parameters.MS.par0(:,k) = par{i}.parameters.MS.par0;
        parameters.MS.logPost0(k) = par{i}.parameters.MS.logPost0;
        parameters.MS.logPost(k) = par{i}.parameters.MS.logPost;
        parameters.MS.gradient(:,k) = par{i}.parameters.MS.gradient;
        parameters.MS.hessian(:,:,k) = par{i}.parameters.MS.hessian;
        parameters.MS.n_objfun(k) = par{i}.parameters.MS.n_objfun;
        parameters.MS.n_iter(k) = par{i}.parameters.MS.n_iter;
        parameters.MS.t_cpu(k) = par{i}.parameters.MS.t_cpu;
        parameters.MS.exitflag(k) = par{i}.parameters.MS.exitflag;
        k=k+1;
    catch
        % catch starts that were not saved
        parameters.MS.par(:,k) = nan(size(parameters.MS.par(:,1)));
        parameters.MS.par0(:,k) = nan(size(parameters.MS.par0(:,1)));
        parameters.MS.logPost0(k) = nan;
        parameters.MS.logPost(k) = nan;
        parameters.MS.gradient(:,k) = nan(size(parameters.MS.gradient(:,1)));
        parameters.MS.hessian(:,:,k) = nan(size(parameters.MS.hessian(:,:,1)));
        parameters.MS.n_objfun(k) = nan;
        parameters.MS.n_iter(k) = nan;
        parameters.MS.t_cpu(k) = nan;
        parameters.MS.exitflag(k) = nan;
        k=k+1;
    end
end

save('parametersToy','parameters','options')