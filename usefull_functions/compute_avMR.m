%% simple average, distance, all sections
ff = gsw_f(nanmean(lt_front));
nobs = nansum(~isnan(glid_MR.CT_io+glid_MR.SA_io),3);
for varr = {'CT','SA','VG','SIG'}
eval([varr{1} '_mean=nanmean(glid_MR.' varr{1} '_io,3);d' varr{1} '_mean=nanstd(glid_MR.' varr{1} '_io,1,3);' varr{1} '_err=d' varr{1} '_mean./sqrt(nobs);'])
end


%% simple average z,distance at seasonal scale
month_front = str2num(datestr(tt_front,'mm'));
iqJFMA = find(month_front==1 | month_front==2 | month_front==3 | month_front==4);
iqMJJA = find(month_front==5 | month_front==6 | month_front==7 | month_front==8);
iqSOND = find(month_front==9 | month_front==10 | month_front==11 | month_front==12);
for varr = {'CT','SA','VG','SIG'}
k=0;
for seas = {'JFMA','MJJA','SOND'}
k=k+1;
eval(['iq_seas = iq' seas{1} ';'])
nobs_seas(:,:,k) = nansum(~isnan(glid_MR.CT_io(:,:,iq_seas)+glid_MR.SA_io(:,:,iq_seas)),3);
eval([lower(varr{1}) '_seas(:,:,k)=nanmean(glid_MR.' varr{1} '_io(:,:,iq_seas),3);d' lower(varr{1}) '_seas(:,:,k)=nanstd(glid_MR.' varr{1} '_io(:,:,iq_seas),1,3);' lower(varr{1}) '_err_seas(:,:,k)=d' lower(varr{1}) '_seas(:,:,k)./sqrt(nobs_seas(:,:,k));'])
end
iqbad = find(nobs_seas<4);
eval([lower(varr{1}) '_seas(iqbad)=NaN;d' lower(varr{1}) '_seas(iqbad)=NaN;' lower(varr{1}) '_err_seas(iqbad)=NaN;'])
end

%% average per season, then annual
for varr = {'CT','SA','VG','SIG'}
eval([lower(varr{1}) '_mean=mean(' lower(varr{1}) '_seas,3);d' lower(varr{1}) '_mean=mean(d' lower(varr{1}) '_seas,3);' lower(varr{1}) '_err=mean(d' lower(varr{1}) '_seas,3);'])
end
