function [mt] = mean_av(t,index)

for i = 1:length(t)
	if ~isnan(t(i))
        mt(i) = nanmean(t(max([1 round(i-index/2)]):min([length(t) round(i+index/2)])));
	else
	mt(i) = NaN;
	end
end
