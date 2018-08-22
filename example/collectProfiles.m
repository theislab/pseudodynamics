clear;
% put path to results of individual profiles here
path = pwd;
par{1} = load([pwd '/parametersProfileExample_' num2str(1) '.mat']);
parameters = par{1}.parameters;
k=2;
for i=2:27
    try
        par{k} = load([pwd '/parametersProfileExample_' num2str(i) '.mat']);
        parameters.P(k).par = par{k}.parameters.P(k).par;
        parameters.P(k).logPost = par{k}.parameters.P(k).logPost;
        parameters.P(k).R = par{k}.parameters.P(k).R;
        parameters.P(k).t_cpu = par{k}.parameters.P(k).t_cpu;
        parameters.P(k).optSteps = par{k}.parameters.P(k).optSteps;
        parameters.P(k).intSteps = par{k}.parameters.P(k).intSteps;
        parameters.P(k).reOptSteps = par{k}.parameters.P(k).reOptSteps;
        k=k+1;
    catch
        parameters.P(k).par = nan(parameters.number,1);
        parameters.P(k).logPost = nan;
        parameters.P(k).R = nan;
        parameters.P(k).t_cpu = nan;
        parameters.P(k).optSteps = nan;
        parameters.P(k).intSteps = nan;
        parameters.P(k).reOptSteps = nan;
        k=k+1;
    end
end
