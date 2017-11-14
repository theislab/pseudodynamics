function [y_a] = augment_cdf(x,x_a,y)
% augment_cdf augments the cdf values y given on a grid x to a finer grid x_a.
i_shift = 0;
y_a = y;
for i = 1:length(x)-1
    i1 = length(find(x_a > x(i)))+length(find(x_a<x(i+1)))-length(x_a); 
    if i1 > 0
        y_a = [y_a(1:i+i_shift,:); bsxfun(@times,y_a(i+i_shift,:),ones(i1,size(y_a,2)));y_a(i+i_shift+1:end,:)];
        i_shift = i_shift+i1;
    end
end
ia = find(x_a < x(1));
y_a = [zeros(length(ia),size(y_a,2));y_a];
ie = find(x_a>x(end));
y_a= [y_a;ones(length(ie),size(y_a,2))];
end

