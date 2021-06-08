
%%%%% independance of estimates
dt = 5;dx = 25;ind_dep = [];
for ll=1:length(tt_front)
dtfront = abs(tt_front-tt_front(ll));
iqt = find(dtfront<dt & dtfront~=0);
if ~isempty(iqt)
for l=1:length(iqt)
dxfront=sw_dist([lt_front(ll) lt_front(iqt(l))],[lg_front(ll) lg_front(iqt(l))],'km');
if ~isempty(dxfront)
disp(['### found cross-section at <' num2str(dt) 'days and <' num2str(dx) 'km : #' num2str(ll) ' and #' num2str(iqt(l)) ])
ind_dep = [[ll iqt(l)];ind_dep];
end
end
end
end
for l=1:size(ind_dep,1)
try ind_dep(find(ind_dep(:,2)==ind_dep(l,1) & ind_dep(:,1)==ind_dep(l,2)),:)=[]; end
end
for l=1:size(ind_dep,1)
try ind_dep(find(ind_dep(:,1)==ind_dep(l,2)),:)=[]; end
end

seas = {'JFMA','MJJA','SOND'};
%%% vmax
clear vmin vbtmin vbcmin
for ll=1:size(glid_MR.VG_io,3)
Dmin = dd_grid(1,glid_MR.iiq1(ll)); Dmax = dd_grid(1,glid_MR.iiq2(ll)); 
tmpbt = glid_MR.VGbt_io(:,:,ll).*double(dd_grid>Dmin & dd_grid<Dmax);tmpbc = glid_MR.VGbc_io(:,:,ll).*double(dd_grid>Dmin & dd_grid<Dmax);tmp = tmpbc+tmpbt;
[a b] = nanmin(tmp);[velm iq]=nanmin(a);
vmin(ll) = velm;vbtmin(ll) = nanmin(tmpbt(:,iq));vbcmin(ll) = nanmin(tmpbc(:,iq));
disp(['#' num2str(ll) ' max vel : ' num2str(velm,3) 'm/s @ ' num2str(dd_grid(b(iq),iq)) 'km, ' num2str(pp_grid(b(iq),iq)) 'm, geostr shear comp ' num2str(vbcmin(ll),2) 'm/s, 1000m comp ' num2str(vbtmin(ll),2) 'm/s'])
end


iqJFMA = find(month_front==1 | month_front==2 | month_front==3 | month_front==4);
iqMJJA = find(month_front==5 | month_front==6 | month_front==7 | month_front==8);
iqSOND = find(month_front==9 | month_front==10 | month_front==11 | month_front==12);
for ll=1:3
eval(['iq_seas = iq' seas{ll} ';'])
disp([seas{ll} ' : Vmax = ' num2str(nanmean(vmin(iq_seas)),2) '+/-' num2str(nanstd(vmin(iq_seas)),2) 'm/s, geostrophic shear ' num2str(nanmean(vbcmin(iq_seas)),2) '+/-' num2str(nanstd(vbcmin(iq_seas)),2) 'm/s, 1000m comp ' num2str(nanmean(vbtmin(iq_seas)),2) '+/-' num2str(nanstd(vbtmin(iq_seas)),2) 'm/s'])
end

%%%% width
width_int = [[(glid_MR.iiq2-glid_MR.iiq1)*(dd_grid(1,2)-dd_grid(1,1))] ];
fit_gauss = [glid_MR.fit_gauss(:,3)' ];
disp(['width_int = ' num2str(nanmean(width_int),2) '+/-' num2str(nanstd(width_int),2) 'km'])
disp(['fit_gauss sigma = ' num2str(nanmean(fit_gauss),2) '+/-' num2str(nanstd(fit_gauss),2) 'km'])
for ll=1:3
eval(['iq_seas = iq' seas{ll} ';'])
disp([seas{ll} ' : width_int = ' num2str(nanmean(width_int(iq_seas)),2) '+/-' num2str(nanstd(width_int(iq_seas)),2) 'km'])
disp([seas{ll} ' : fit_gauss sigma = ' num2str(nanmean(fit_gauss(iq_seas)),2) '+/-' num2str(nanstd(fit_gauss(iq_seas)),2) 'km'])
end

%%%%
%%%%% synthesis transport
for transp = {'','bt','bc'} % mass heat salt
for varr = {'','H','S'} % mass heat salt
%%% different water masses
for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
eval([transp{1} varr{1} 'Tr' watermass{1} ' = -[TR.' transp{1} varr{1} 'Tr' watermass{1} '];'])
%' transp{1} varr{1} 'Tr' watermass{1} '([20 23]) = -TR.' transp{1} varr{1} 'Tr' watermass{1} '_ideal([20 23]);'])
end
end
end
for varr = {'tt_front','width_int','vmin','vbtmin','vbcmin','fit_gauss','Tr','btTr','bcTr','TrAW','btTrAW','bcTrAW','HTrAW','TrmAW','TrnAW','TrsASW','TrwASW'}
eval(['tmp = ' varr{1} ';'])
for l=1:size(ind_dep,1)
eval(['tmp(ind_dep(l,:))=[];tmp=[tmp nanmean(' varr{1} '(ind_dep(l,:)))];'])
end
eval([varr{1} ' = tmp;'])
end

month_front = str2num(datestr(tt_front,'mm'));
iqJFMA = find(month_front==1 | month_front==2 | month_front==3 | month_front==4);
iqMJJA = find(month_front==5 | month_front==6 | month_front==7 | month_front==8);
iqSOND = find(month_front==9 | month_front==10 | month_front==11 | month_front==12);

seasweight = ones(size(vmin));
for ll=1:3
eval(['seasweight(iq' seas{ll} ') = 1/length(iq' seas{ll} ');'])
end

%%%% make latex table
disp('\begin{table}')
disp('\small')
disp('\begin{tabular}{lc|ccc|c}')
disp(' & units & Winter (JFMA) & Summer (MJJA) & Fall (SOND) & Annual \\ \hline')
disp(['$n$ &  & ' num2str(length(iqJFMA)) ' & ' num2str(length(iqMJJA)) ' & ' num2str(length(iqSOND)) ' & ' num2str(length(month_front)) ' \\'])
%disp(['$\partial_x T$, max  & ' num2str(gradmax,3) ' & ' num2str(gradmax,3) ' & ' num2str(gradmax,3) ' & ' num2str(gradmax,3) ' \\'])
disp(['$U_{max}$ & m\,s$^{-1}$ & ' num2str(nanmean(-vmin(iqJFMA)),'%1.2f') '$\pm$' num2str(nanstderr(vmin(iqJFMA)),'%1.2f') ' & ' num2str(nanmean(-vmin(iqMJJA)),'%1.2f') '$\pm$' num2str(nanstderr(vmin(iqMJJA)),'%1.2f') ' & ' num2str(nanmean(-vmin(iqSOND)),'%1.2f') '$\pm$' num2str(nanstderr(vmin(iqSOND)),'%1.2f') ' & ' num2str(wmean(-vmin,seasweight),'%1.2f') '$\pm$' num2str(nanstderr(vmin),'%1.2f') ' \\'])
disp(['$U^{shear}_{max}$ $\vert$ $U^{1000}_{max}$ & m\,s$^{-1}$ & ' num2str(nanmean(-vbcmin(iqJFMA)),'%1.2f') '$\pm$' num2str(nanstderr(vbcmin(iqJFMA)),'%1.2f') ' $\vert$ ' num2str(nanmean(-vbtmin(iqJFMA)),'%1.2f') '$\pm$' num2str(nanstderr(vbtmin(iqJFMA)),'%1.2f') ' & ' num2str(nanmean(-vbcmin(iqMJJA)),'%1.2f') '$\pm$' num2str(nanstderr(vbcmin(iqMJJA)),'%1.2f') ' $\vert$ ' num2str(nanmean(-vbtmin(iqMJJA)),'%1.2f') '$\pm$' num2str(nanstderr(vbtmin(iqMJJA)),'%1.2f') ' & ' num2str(nanmean(-vbcmin(iqSOND)),'%1.2f') '$\pm$' num2str(nanstderr(vbcmin(iqSOND)),'%1.2f') ' $\vert$ ' num2str(nanmean(-vbtmin(iqSOND)),'%1.2f') '$\pm$' num2str(nanstderr(vbtmin(iqSOND)),'%1.2f') ' & ' num2str(wmean(-vbcmin,seasweight),'%1.2f') '$\pm$' num2str(nanstderr(vbcmin),'%1.2f') ' $\vert$ ' num2str(wmean(-vbtmin,seasweight),'%1.2f') '$\pm$' num2str(nanstderr(vbtmin),'%1.2f') '\\'])
disp(['$\sigma_{fit}$ & km & ' num2str(nanmean(fit_gauss(iqJFMA)),'%2.0f') '$\pm$' num2str(nanstderr(fit_gauss(iqJFMA)),'%2.0f') ' & ' num2str(nanmean(fit_gauss(iqMJJA)),'%2.0f') '$\pm$' num2str(nanstderr(fit_gauss(iqMJJA)),'%2.0f') ' & ' num2str(nanmean(fit_gauss(iqSOND)),'%2.0f') '$\pm$' num2str(nanstderr(fit_gauss(iqSOND)),'%2.0f') ' & ' num2str(wmean(fit_gauss,seasweight),'%2.0f') '$\pm$' num2str(nanstderr(fit_gauss),'%2.0f') ' \\'])
disp(['$L_{int}$ & km & ' num2str(nanmean(width_int(iqJFMA)),'%2.0f') '$\pm$' num2str(nanstderr(width_int(iqJFMA)),'%2.0f') ' & ' num2str(nanmean(width_int(iqMJJA)),'%2.0f') '$\pm$' num2str(nanstderr(width_int(iqMJJA)),'%2.0f') ' & ' num2str(nanmean(width_int(iqSOND)),'%2.0f') '$\pm$' num2str(nanstderr(width_int(iqSOND)),'%2.0f') ' & ' num2str(wmean(width_int,seasweight),'%2.0f') '$\pm$' num2str(nanstderr(width_int),'%2.0f') ' \\'])
disp(['$Tr_{AW}$ & Sv & ' num2str(nanmean(TrAW(iqJFMA)),'%1.1f') '$\pm$' num2str(nanstderr(TrAW(iqJFMA)),'%1.1f') ' & ' num2str(nanmean(TrAW(iqMJJA)),'%1.1f') '$\pm$' num2str(nanstderr(TrAW(iqMJJA)),'%1.1f') ' & ' num2str(nanmean(TrAW(iqSOND)),'%1.1f') '$\pm$' num2str(nanstderr(TrAW(iqSOND)),'%1.1f') ' & ' num2str(wmean(TrAW,seasweight),'%1.1f') '$\pm$' num2str(nanstderr(TrAW),'%1.1f') ' \\'])
disp(['$Tr^{shear}_{AW}$ $\vert$ $Tr^{1000}_{AW}$ & Sv & ' num2str(nanmean(bcTrAW(iqJFMA)),'%1.1f') '$\pm$' num2str(nanstderr(bcTrAW(iqJFMA)),'%1.1f') ' $\vert$ ' num2str(nanmean(btTrAW(iqJFMA)),'%1.1f') '$\pm$' num2str(nanstderr(btTrAW(iqJFMA)),'%1.1f') ' & ' num2str(nanmean(bcTrAW(iqMJJA)),'%1.1f') '$\pm$' num2str(nanstderr(bcTrAW(iqMJJA)),'%1.1f') ' $\vert$ ' num2str(nanmean(btTrAW(iqMJJA)),'%1.1f') '$\pm$' num2str(nanstderr(btTrAW(iqMJJA)),'%1.1f') ' & ' num2str(nanmean(bcTrAW(iqSOND)),'%1.1f') '$\pm$' num2str(nanstderr(bcTrAW(iqSOND)),'%1.1f') ' $\vert$ ' num2str(nanmean(btTrAW(iqSOND)),'%1.1f') '$\pm$' num2str(nanstderr(btTrAW(iqSOND)),'%1.1f') ' & ' num2str(wmean(bcTrAW,seasweight),'%1.1f') '$\pm$' num2str(nanstderr(bcTrAW),'%1.1f') ' $\vert$ ' num2str(wmean(btTrAW,seasweight),'%1.1f') '$\pm$' num2str(nanstderr(btTrAW),'%1.1f') ' \\'])
disp(['$HTr_{AW}$ & TW & ' num2str(round(nanmean(HTrAW(iqJFMA))/1e12),'%3.0f') '$\pm$' num2str(round(nanstderr(HTrAW(iqJFMA))/1e12),'%3.0f') ' & ' num2str(round(nanmean(HTrAW(iqMJJA))/1e12),'%3.0f') '$\pm$' num2str(round(nanstderr(HTrAW(iqMJJA))/1e12),'%3.0f') ' & ' num2str(round(nanmean(HTrAW(iqSOND))/1e12),'%3.0f') '$\pm$' num2str(round(nanstderr(HTrAW(iqSOND))/1e12),'%3.0f') ' & ' num2str(round(wmean(HTrAW,seasweight)/1e12),'%3.0f') '$\pm$' num2str(round(nanstderr(HTrAW)/1e12),'%3.0f') ' \\'])
disp('\end{tabular}')
disp('\caption{Seasonal and mean Polar Front characteristics (average $\pm$ standard error) : number of independent sections ($n$), maximum absolute frontal velocity ($U_{max}$) decomposed into maximum velocity due to the geostrophic shear above 1000\,m ($U^{shear}_{max}$) and maximum velocity at 1000\,m ($U^{1000}_{max}$), half-width of the fitted Gaussian function ($\sigma_{fit}$), width of the integration window ($L_{int}$), AW transport ($Tr_{AW}$) decomposed into AW transport due to geostrophic shear above 1000\,m ($Tr^{shear}_{AW}$) and AW transport due to absolute velocity at 1000\,m ($Tr^{1000}_{AW}$), AW heat transport ($HTr_{AW}$). For the annual mean, a weighted average according to the number of sections performed during each season is considered.}')
disp('\label{tab-front}')
disp('\end{table}')
