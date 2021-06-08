function cost = coast_gaussian_MRjet(param)
global Djet Vjet iqtest

Vtest = nanmin(Vjet(iqtest))*gauss_distribution(Djet, param(1), param(3));

cost = nanmean( (Vjet(iqtest) - Vtest(iqtest)).^2 );
