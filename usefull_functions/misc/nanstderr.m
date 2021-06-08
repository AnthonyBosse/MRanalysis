function [y] = nanstderr(x)
y = nanstd(x)./sqrt(length(find(~isnan(x))));
