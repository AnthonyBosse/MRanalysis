function [iq_sASW iq_wASW iq_AW iq_mAW iq_nAW iq_NSIW iq_NSDW] = select_water_mass(sa,ct)

sig0 = gsw_sigma0(sa,ct);

iq_wASW = (sig0>27.8 & sa<35.07);iq_sASW = (sa<35.12 & sig0<27.8);
iq_AW = (sa>35.17);iq_nAW = (sa>35.27 & ct>6);iq_mAW = (sig0<27.82 & sig0>27.71 & sa>35.28 & sa<35.34);
iq_NSIW = (sig0>28 & sa<35.08 & sa>35.07);iq_NSDW = (sig0>28 & sa>35.08 & sa<35.09);

%%%%%
if 0
col_paired = cptcmap('GMT_paired');
ctdt_sum = cruise_MR(1).CTD_CT(:,ind_summer);ctds_sum = cruise_MR(1).CTD_SA(:,ind_summer);
ctdt_wint = cruise_MR(2).CTD_CT(:,ind_winter);ctds_wint = cruise_MR(2).CTD_SA(:,ind_winter);
iqMR = find(lg_all<4.5 & lt_all>70.75);
glidt = CT_all(:,iqMR);glids = SA_all(:,iqMR);

[mz,xmid,ymid,numz,stdz]=twodstats([ctds_sum(:)' ctds_wint(:)' glids(:)'],[ctdt_sum(:)' ctdt_wint(:)' glidt(:)'],ones(size([ctds_sum(:)' ctds_wint(:)' glids(:)'])),200);
numz(numz==0)=NaN;
[sgrid tgrid] = meshgrid(xmid,ymid);
sig0 = gsw_sigma0(sgrid,tgrid);sig1 = gsw_sigma1(sgrid,tgrid);

close all
axes('position',[.1 .1 .5 .6])
hold on;grid on;box on
pcolor(sgrid,tgrid,log10(numz));shading flat
[c,h] = contour(sgrid,tgrid,sig0,[26:.1:28.15],'k')
clabel(c,h,'fontsize',6);
caxis([-1 4])
cptcmap('stern','flip',true)
h=colorbarnew('v',0.02,1,'# count');set(h,'ylim',[0 4]);set(h,'ytick',[0:4],'yticklabel',{'1';'10';'10^2';'10^3';'10^4'});
xlabel('S_A (g kg^{-1})');ylabel('\Theta (\circC)')
%% water mass
i_wASW = sig0>27.8 & sgrid<35.07;contour(sgrid,tgrid,double(i_wASW),[0.5 0.5],'color',col_paired(2,:),'linewidth',1)
i_sASW = sgrid<35.12 & sig0<27.8;contour(sgrid,tgrid,double(i_sASW),[0.5 0.5],'color',col_paired(2,:),'linewidth',1)

i_AW = sgrid<35.17;contour(sgrid,tgrid,double(i_AW),[0.5 0.5],'-','linewidth',1,'color',col_paired(5,:))
i_nAW = sgrid>35.27 & tgrid>6;contour(sgrid,tgrid,double(i_nAW),[0.5 0.5],'--','color',col_paired(6,:),'linewidth',1)
i_mAW = sig0<27.82 & sig0>27.71 & sgrid>35.28 & sgrid<35.34;contour(sgrid,tgrid,double(i_mAW),[0.5 0.5],'-','color',col_paired(7,:),'linewidth',1)

i_NSDW = sig0>28 & sgrid>35.08 & sgrid<35.09;contour(sgrid,tgrid,double(i_NSDW),[0.5 0.5],'color',col_paired(10,:),'linewidth',1)
i_NSIW = sig0>28 & sgrid<35.08 & sgrid>35.07;contour(sgrid,tgrid,double(i_NSIW),[0.5 0.5],'color',col_paired(9,:),'linewidth',1)
%% labels
text(35.18,10,'AW','color',col_paired(5,:),'fontsize',9);text(35.33,10,'nAW','color',col_paired(6,:),'fontsize',9);text(35.3,4.5,'mAW','color',col_paired(7,:),'fontsize',9);
text(35.09,-0.6,'NSDW','color',col_paired(10,:),'fontsize',9);text(35.01,-0.6,'NSIW','color',col_paired(9,:),'fontsize',9);
text(35,1.5,'wASW','color',col_paired(2,:),'fontsize',9);text(35.05,7.5,'sASW','color',col_paired(2,:),'fontsize',9);
%%
printHR('TS_watermass')
end
