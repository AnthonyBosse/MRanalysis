function [mld_tp pref] = compute_mld(p,tp,delta_tp)

%if nanmax(tp(:))-nanmin(tp(:))<0.05
%delta_tp = 0.01;
%end

mld_tp = NaN;
iqnan = find(~isnan(p+tp));
[ps I J] = unique(p(iqnan));
if length(ps)>1

if nanmean(diff(ps))>5
pmin = round(nanmin(p(find(~isnan(tp)))));
tpp_interp = mean_av(interp1(ps,tp(iqnan(I)),pmin:2700),10);pp_int = pmin:2700;
else
tpp_interp = tp(iqnan(I));pp_int = ps;
end

if length(pp_int)>100/nanmean(diff(ps))

% ref at 10m (ou au delà)
iq = find(pp_int<=10);
if isempty(iq)
iq = nanmin(find(~isnan(tpp_interp)));
if pp_int(iq)>50
iq = [];
end
end

pref = nanmean(pp_int(iq));
if ~isempty(iq)
%X = robust_median(tpp_interp(iq)+rand(size(iq))*0.0000001,1);

% old method wrong!!
%tp_ref = nanmean(tpp_interp(iq));
%mld_tp = 10+pp_int(nanmin(find(abs(tpp_interp(10:end)-tp_ref)>delta_tp)));


tp_ref = nanmean(tpp_interp(iq));
iq_out = 1:length(pp_int);iq_out = setdiff(iq_out,iq);
mld_tp = pp_int(nanmin(find(abs(tpp_interp-tp_ref)>delta_tp)));

%if ~isempty(mld_tp) & mld_tp>300 & length(tpp_interp)>305
%tp_ref = nanmean(tpp_interp(295:305));
%tpp_interp(pp_int<300) = NaN;
%mld_tp = pp_int(nanmin(find(abs(tpp_interp-tp_ref)>delta_tp/3)));
%end
if isempty(mld_tp)
mld_tp = nanmax(pp_int(find(~isnan(tpp_interp))));
end

end
end
end

if isnan(mld_tp); pref = NaN; end

check_plot = 0;
if check_plot
figure
hold on;grid on;box on
plot(tp,p,'.','markersize',4)
plot(tpp_interp,pp_int,'-')
plot([tp_ref tp_ref],[0 15],'-or')
plot([tp_ref+delta_tp tp_ref+delta_tp],[0 nanmax(p)],'--r')
plot(interp1(pp_int,tpp_interp,mld_tp),mld_tp,'vg')
zdn
end
