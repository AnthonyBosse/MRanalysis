% load ERA5 24h and get wind stress curl over LB
load ERA5_PROVOLO_24h_tau.mat
bat_era5 = interp2(lg_bat,lt_bat,z_smo,lon,lat);
mask_LB = double(lon>-1 & lon<20 & lat>68.25 & lat<72.25 & bat_era5<-2800);mask_LB(mask_LB==0)=NaN;
[xera yera] = earthcoord('flat',lon,lat,nanmean(lon(:)),nanmean(lat(:)));
for l=1:size(u10_av,1)
curl_tau(l,:,:)=flipud(spherecurl(flipud(lat(:,1)),lon(1,:),flipud(squeeze(taux_av(l,:,:))),flipud(squeeze(tauy_av(l,:,:))))); 
tmp = squeeze(curl_tau(l,:,:)).*mask_LB; mean_curltau_LB(l) = nanmean(tmp(:));
end
surf_LB=(xera(1,2)-xera(1,1))*(yera(1,1)-yera(2,1))*1e6*nansum(mask_LB(:)); % surface LB in m2
[a b]=contour(xera*1000,yera*1000,double(~isnan(mask_LB)),[.5 .5]);
xLB = a(1,2:end-6); yLB = a(2,2:end-6);
Req_LB = sqrt(surf_LB/pi);
circ_LB = nansum(sqrt(diff(xLB).^2+diff(yLB).^2));


ttmonthly = [datenum(2016,6:12,1) datenum(2017,1:8,1)];
bincurl = bin1d(tt_av,mean_av(mean_curltau_LB,45)',j*ttmonthly);
binTrAW = bin1d(tt_front,TrAW,j*ttmonthly);
binbtTrAW = bin1d(tt_front,btTrAW,j*ttmonthly);

%%% figure tseries
lab_l = {'\bfa','\bfb','\bfc','\bfd','\bfe','\bff','\bfg','\bfh','\bfi','\bfj','\bfk','\bfl'};
ylab_l = {'Tr (Sv)','Tr (Sv)','\nabla\times\tau (10^{-7} N m^{-3})','V_{max} (m s^{-1})'};
seas_l = {'JFMA','MJJA','SOND'};
tit_l = {'\Theta (\circC)','S_A (g kg^{-1})','V_g (m s^{-1})','Obs. Probability'};
ncol = 1;nlig = 3;
close all
[ax,row,col,IsTop,IsBottom,IsLeft,IsRight]=subaxes(nlig,ncol,'RowSpacing',0.2,'ColSpacing',0.1,'PageMargins',[0.08 0.28 0.07 0.01]);

for l=1:length(ax)
if ~IsLeft(l)
set(ax(l),'yticklabel','')
else
ylabel(ax(l),ylab_l{row(l)})
end
if ~IsBottom(l)
set(ax(l),'xticklabel','')
end
hold(ax(l),'on');grid(ax(l),'on');box(ax(l),'on');
end


l=1;TrASW = TrwASW+TrsASW;TrDW = TrNSIW+TrNSDW;
a=bar(ax(l),tt_front,TrAW,1.5);set(a,'edgecolor',[1 1 1],'facecolor',col_paired(5,:));
c=bar(ax(l),tt_front,TrmAW,1.5);set(c,'edgecolor',[1 1 1],'facecolor',col_paired(8,:));
b=bar(ax(l),tt_front,TrnAW,0.5);set(b,'edgecolor','none','facecolor',col_paired(6,:));
d=bar(ax(l),tt_front,TrASW,1);set(d,'edgecolor','none','facecolor',col_paired(2,:));
set(ax(l),'ylim',[0 nanmax(TrAW).*1.1],'xlim',mima(ttmonthly),'xtick',ttmonthly)
leg=legend([a b c d],{'AW','nAW','mAW','ASW'},'numcolumns',1,'location','northwest','box','off');pos=get(leg,'position');set(leg,'position',[pos(1)+0.01 pos(2)+0.04 pos(3) pos(4)*.5])
labfig('\bfa',1,0,0,ax(l));


l=2;
a=bar(ax(l),tt_front,TrAW,1.5);set(a,'edgecolor',[1 1 1],'facecolor',col_paired(5,:));
b=bar(ax(l),tt_front,bcTrAW,1.5);set(b,'edgecolor',[1 1 1],'facecolor',[.5 .5 .5],'facealpha',.5);
c=bar(ax(l),tt_front,TrAW-bcTrAW,0.5);set(c,'edgecolor','none','facecolor',[0 0 0]);
set(ax(l),'ylim',[0 nanmax(TrAW).*1.1],'xlim',mima(ttmonthly),'xtick',ttmonthly)
leg=legend([a b c],{'total','geostr. shear','deep curr.'},'numcolumns',1,'location','northwest','box','off');pos=get(leg,'position');set(leg,'position',[pos(1)+0.01 pos(2)+0.035 pos(3) pos(4)*.5])
labfig('\bfb',1,0,0,ax(l));

l=3;
plot(ax(l),ttmonthly,zeros(size(ttmonthly)),'-k')
[htop,hbot] = anomaly(tt_av,mean_av(mean_curltau_LB*1e7,30),'topcolor','r','bottomcolor','b','edgecolor','none');
set(htop,'parent',ax(l));set(hbot,'parent',ax(l));
st=stairs(ax(l),[ttmonthly],[bincurl.mean bincurl.mean(end)]*1e7,'-','linewidth',1,'color','k');
set(ax(l),'ylim',mima(mean_av(mean_curltau_LB*1e7,30)).*1.1,'xlim',mima(ttmonthly),'xtick',ttmonthly)
dtickx
labfig('\bfc',1,0,0,ax(l));

pos = get(ax(l),'position');axes('position',[pos(1) pos(2)-0.05 pos(3:4)])
hold on; axis off
text(datenum(2016,9,1),0,'2016','fontsize',6)
text(datenum(2017,1,1),0,'|','fontsize',6)
text(datenum(2017,4,1),0,'2017','fontsize',6)
xlim([datenum(2016,6,1) datenum(2017,8,1)])

%% seasonal cycle
[ax2,row2,col2,IsTop2,IsBottom2,IsLeft2,IsRight2]=subaxes(nlig,ncol,'RowSpacing',0.2,'ColSpacing',0.1,'PageMargins',[0.82 0.01 0.07 0.01]);
for l=2:length(ax)
if ~IsLeft(l)
set(ax2(l),'yticklabel','')
else
ylabel(ax2(l),ylab_l{row(l)})
end
if ~IsBottom(l)
set(ax2(l),'xticklabel','')
end
hold(ax2(l),'on');grid(ax2(l),'on');box(ax2(l),'on');
set(ax2(l),'xlim',[0.5 12.5],'xtick',[1:1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'})
end
axis(ax2(1),'off');


%% fit seasonal %%% fit : B sin(2*pi*(x-C)/12)
B=0:.1:3;C=1:1:12;
for l=1:length(B)
for k=1:length(C)
cost_tot(l,k) = nansum( (TrAW-nanmean(TrAW) - B(l)*cos(2*pi*(month_front-0.5-C(k))/12)').^2 );
cost_bc(l,k) = nansum( (bcTrAW-nanmean(bcTrAW) - B(l)*cos(2*pi*(month_front-0.5-C(k))/12)').^2 );
cost_bt(l,k) = nansum( (TrAW-bcTrAW-nanmean(TrAW-bcTrAW) - B(l)*cos(2*pi*(month_front-0.5-C(k))/12)').^2 );
end
end
l=2;
[a b] = nanmin(cost_tot);[c d] = nanmin(a); Btot = B(b(d));Ctot = C(d); TrAW_seasfit = nanmean(TrAW) + Btot*cos(2*pi*((1:12)-0.5-Ctot)/12);
[a b] = nanmin(cost_bc);[c d] = nanmin(a); Bbc = B(b(d));Cbc = C(d); TrbcAW_seasfit = nanmean(bcTrAW) + Bbc*cos(2*pi*((1:12)-0.5-Cbc)/12);
[a b] = nanmin(cost_bt);[c d] = nanmin(a); Bbt = B(b(d));Cbt = C(d); TrbtAW_seasfit = nanmean(TrAW-bcTrAW) + Bbt*cos(2*pi*((1:12)-0.5-Cbt)/12);
plot(ax2(l),1:12,TrAW_seasfit,'-','color',col_paired(5,:),'linewidth',1);
plot(ax2(l),1:12,TrbcAW_seasfit,'-','color',nanmean([.5 .5 .5;col_paired(5,:)]),'linewidth',1);
plot(ax2(l),1:12,TrbtAW_seasfit,'-','color','k','linewidth',1);

a=plot(ax2(l),month_front,TrAW,'ow','markerfacecolor',col_paired(5,:),'markersize',4);
a=plot(ax2(l),month_front,bcTrAW,'ow','markerfacecolor',nanmean([.5 .5 .5;col_paired(5,:)]),'markersize',4);
a=plot(ax2(l),month_front,TrAW-bcTrAW,'ok','markersize',2,'markerfacecolor','k');
set(ax2(l),'ylim',[0 nanmax(TrAW)].*1.1)
labfig('\bfd',3,0,0,ax2(l));


l=3;
c=bar(ax2(l),1:12,mean_curltau_LB_monthly*1e7,1);set(c,'edgecolor','w','facecolor',[.5 .5 .5]);
iq=find(str2num(datestr(bincurl.x,'yyyy'))==2016);
aa=plot(ax2(l),str2num(datestr(bincurl.x(iq),'mm')),bincurl.mean(iq)*1e7,'ow','markerfacecolor',get(st,'color'));
iq=find(str2num(datestr(bincurl.x,'yyyy'))==2017);
bb=plot(ax2(l),str2num(datestr(bincurl.x(iq),'mm')),bincurl.mean(iq)*1e7,'sw','markerfacecolor',get(st,'color'));
labfig('\bfe',3,0,0,ax2(l));


printHR('fig/tseriesMR_n')