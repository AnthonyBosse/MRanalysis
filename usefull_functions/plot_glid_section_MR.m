lx = .37;ly = 0.5*lx/21*29.7;llx = 0.08;

[glid_MR.maskio_sASW glid_MR.maskio_wASW glid_MR.maskio_AW glid_MR.maskio_mAW glid_MR.maskio_nAW glid_MR.maskio_NSIW glid_MR.maskio_NSDW] = select_water_mass(glid_MR.SA_io,glid_MR.CT_io);
glid_MR.maskio_ASW = ceil((glid_MR.maskio_wASW+glid_MR.maskio_sASW)/2);

%%%%%% spatial, d cross front,
m_proj('equidist','lat',lat_MR,'lon',lon_MR);

close all;figure('visible',opt_vis)
axes('position',[llx-0.02 .55 lx lx],'xtick',[-200:50:200],'xlim',[-100 100],'ylim',[-20 5])
hold on;grid on;box on
m_grid('linestyle','none','tickdir','out')
m_contour(lg_bat,lt_bat,z,[-3500:100:500],'color',[.8 .8 .8])
m_contour(lg_bat,lt_bat,z,[-3500:500:500],'color',[.2 .2 .2])
m_plot(lg_all,lt_all,'color','b')
m_scatter(lg_all(ind_front),lt_all(ind_front),20,nanmean(CT_all(300:400,ind_front)),'filled')
indv = find(TVEL_all>nanmin(tt_all(ind_front)) & TVEL_all<nanmax(tt_all(ind_front)));
m_quiver(LGV_all(indv),LTV_all(indv),VX_all(indv),VY_all(indv),2,'r','ShowArrowHead','off')
m_plot(lg_front(ll),lt_front(ll),'ok','markerfacecolor','w','markersize',5)

caxis(t_lim); colormap(cmocean('thermal'))
colorbarnew('v',-.01,1,'\Theta_{300-400} (\circC)');

m_proj('equidist','lat',mima(llt_glid)+[-.2 .2],'lon',mima(llg_glid) + [-.5 .5]);
axes('position',[llx+.45 .55 lx lx],'xtick',[-200:50:200],'xlim',[-100 100],'ylim',[-20 5])
hold on;grid on;box on
m_grid('linestyle','none','tickdir','out')
m_contour(lg_bat,lt_bat,z,[-3500:100:500],'color',[.8 .8 .8])
m_contour(lg_bat,lt_bat,z,[-3500:500:500],'color',[.2 .2 .2])
m_plot(lg_all,lt_all,'color','b')
m_scatter(lg_all(ind_front),lt_all(ind_front),20,nanmean(CT_all(300:400,ind_front)),'filled')
indv = find(TVEL_all>nanmin(tt_all(ind_front)) & TVEL_all<nanmax(tt_all(ind_front)));
m_quiver(LGV_all(indv),LTV_all(indv),VX_all(indv),VY_all(indv),2,'r','ShowArrowHead','off')
m_plot(lg_front(ll),lt_front(ll),'ok','markerfacecolor','w','markersize',5)

caxis(t_lim); colormap(cmocean('thermal'))

axes('position',[llx llx .8 .4])

yyaxis left
plot(tt_all(ind_front),d_sec)
hold on;grid on;box on
plot(tt_all(ind_front),dd_sec)
dtickx
ylabel('d_{front}')

yyaxis right
plot(tt_all(ind_front),nanmean(glid_pv.sint(:,ind_front)))
ylabel('sin(\theta)')

eval(['printHR(''fig/traj_MR' num2str(ll) ''')'])


%%%%%

close all;figure('visible',opt_vis)

ax1=axes('position',[llx llx+ly+0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
pcolor(dd_sec,0:1010,ct_tmp);shading interp
[c,h]=contour(dd_sec,0:1010,sig_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(ax1,cmocean('thermal'));caxis(t_lim);zdn;ylim([0 1000])
ylabel('Depth (m)');
xlabel('Cross-front distance (km)')

ax2=axes('position',[llx+lx+llx-0.01 llx+ly+0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
pcolor(glid_MR.bin_d,glid_MR.bin_p,glid_MR.CT_io(:,:,ll));shading interp
[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(ax2,cmocean('thermal'));caxis(t_lim);zdn;ylim([0 1000])
colorbarnew('v',0.01,1,'\Theta (\circC)');
xlabel('Cross-front distance (km)')


axes('position',[llx llx+2*ly+2*0.02 lx ly],'yticklabel',[],'xtick',[max(nanmin(dd_sec),-150) min(nanmax(dd_sec),150)],'xlim',[-100 100],'xticklabel',{[datestr(tt_all(ind_front([1 end])),'dd mmmm, YYYY')]},'xaxislocation','top');
ax3=axes('position',[llx llx+2*ly+2*0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xticklabel','','xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
pcolor(dd_sec,0:1010,sa_tmp);shading interp
[c,h]=contour(dd_sec,0:1010,sig_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(ax3,cmocean('haline'));caxis(s_lim);zdn;ylim([0 1000])
ylabel('Depth (m)')

ax4=axes('position',[llx+lx+llx-0.01 llx+2*ly+2*0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xticklabel','','xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
pcolor(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SA_io(:,:,ll));shading interp
[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(ax4,cmocean('haline'));caxis(s_lim);zdn;ylim([0 1000])
colorbarnew('v',0.01,1,'S_A (g kg^{-1})');
title('Optimal interpolation')


eval(['printHR(''fig/radTS_MR' num2str(ll) ''')'])



%%%%%%%%%%%%%%%%%%%% Vcyclo, N2, PV
lx = .35;ly = 0.5*lx/21*29.7;llx = 0.12;
v_lim = [-.5 .5];vort_lim = [-0.15 0.15];

close all;figure('visible',opt_vis)


axes('position',[0.08 0.08+ly+0.02 lx ly],'yticklabel',[],'xtick',[max(nanmin(dd_sec),-150) min(nanmax(dd_sec),150)],'xlim',[-100 100],'xticklabel',{[datestr(tt_all(ind_front([1 end])),'dd mmmm, YYYY')]},'xaxislocation','top')
axes('position',[0.08 0.08+ly+0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;grid on;box on
pcolor(dd_sec,glid_pv.pp(:,1),v_tot_tmp);shading interp
[c,h]=contour(dd_sec,glid_pv.pp(:,1),sigs_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
colormap(cmocean('delta'));caxis(v_lim)
zdn;ylim([0 1000]);ylabel('Depth (m)')
xlabel('Cross-front distance (km)')


axes('position',[0.08+lx+llx 0.08+ly+0.02 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;box on;grid on
pcolor(glid_MR.bin_d,glid_MR.bin_p,glid_MR.VG_io(:,:,ll));shading interp
[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
plot([Djet(iiq1) Djet(iiq1) Djet(iiq2) Djet(iiq2) Djet(iiq1)],[0 1000 1000 0 0],'-m');
colormap(cmocean('delta'));caxis(v_lim);
zdn;ylim([0 1000]);ylabel('Depth (m)')
colorbarnew('v',0.02,1,'V_g (m s^{-1})');
title('Optimal interpolation')
xlabel('Cross-front distance (km)')


axes('position',[0.08+lx+llx 0.08+2*ly+0.15 lx ly],'xtick',[-200:50:200],'xlim',[-100 100])
yyaxis right
hold on;grid on;box on
if 0
plot(glid_MR.bin_d(1,:),glid_MR.Tr(ll,:),'-k'); plot(glid_MR.bin_d(1,:),glid_MR.TrAW(ll,:),'-r');
plot(glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll})),glid_MR.Tr(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll})),'-k','linewidth',2);
plot(glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll})),glid_MR.TrAW(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll})),'-r','linewidth',2)
plot([glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))],[glid_MR.Tr(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))],'ow','markerfacecolor','k','markersize',5)
text([glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))]+10,[glid_MR.Tr(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))]-0.5,['\bfTot: ' num2str(round(glid_MR.Tr_tot(ll)*10)/10) ' Sv'],'color','k','fontsize',7)
plot([glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))],[glid_MR.TrAW(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))],'ow','markerfacecolor','r','markersize',5)
text([glid_MR.bin_d(1,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))]+10,[glid_MR.TrAW(ll,glid_MR.indneg{ll}(glid_MR.front_ind{ll}(end)))]+0.5,['\bfAW: ' num2str(round(glid_MR.Tr_AW(ll)*10)/10) ' Sv'],'color','r','fontsize',7)
end

plot(glid_MR.bin_d(1,:),TR.Tr_int_gauss(ll,:),'-k');aa= plot(glid_MR.bin_d(1,:),TR.TrAW_int_gauss(ll,:),'-');
plot(glid_MR.bin_d(1,[iiq1:iiq2]),TR.Tr_int_gauss(ll,[iiq1:iiq2]),'-k','linewidth',2);
plot(glid_MR.bin_d(1,[iiq1:iiq2]),TR.TrAW_int_gauss(ll,[iiq1:iiq2]),'-','linewidth',2,'color',get(aa,'color'))
plot([glid_MR.bin_d(1,iiq2)],[TR.Tr_int_gauss(ll,iiq2)],'ow','markerfacecolor','k','markersize',5)
text([glid_MR.bin_d(1,iiq2)]+10,[TR.Tr_int_gauss(ll,iiq2)]-0.5,['\bfTot: ' num2str(round(TR.Tr(ll)*10)/10) ' Sv'],'color','k','fontsize',7)
plot([glid_MR.bin_d(1,iiq2)],[TR.TrAW_int_gauss(ll,iiq2)],'ow','markerfacecolor',get(aa,'color'),'markersize',5)
text([glid_MR.bin_d(1,iiq2)]+10,[TR.TrAW_int_gauss(ll,iiq2)]+0.5,['\bfAW: ' num2str(round(TR.TrAW(ll)*10)/10) ' Sv'],'color',get(aa,'color'),'fontsize',7)

ylim([-15 2.5]);ylabel('Tr [Sv]');

yyaxis left
aa=plot(glid_MR.bin_d(1,:),vsurf_io);
plot(glid_MR.bin_d(1,indneg_io(front_ind_io)),vsurf_io(indneg_io(front_ind_io)),'-','linewidth',2,'color',get(aa,'color'))
plot(Djet,Vfitg,'--g')
plot([Djet(iiq1) Djet(iiq1)],[-1 1],'--g');plot([Djet(iiq2) Djet(iiq2)],[-1 1],'--g')
ylabel('V_g^{0-25} (m s^{-1})')
ylim([-0.6 0.1])
xlabel('Cross-front distance (km)')

eval(['printHR(''fig/radDiag_MR' num2str(ll) ''')'])



lx = .38;ly = 0.5*lx/21*29.7;llx = 0.08;
%%%%  TS diagram and transport in TS class
[X,Y] = meshgrid(linspace(34.9,35.4,100),linspace(-1,10,100));
Z0 = gsw_sigma0(X,Y);

close all;figure('visible',opt_vis)
axes('position',[llx-0.025 .55 lx*0.9 lx*1.1])
hold on;box on
[c,h]=contour(X,Y,Z0,[26.8:.025:28.4],'color',[.8 .8 .8]);
[c,h]=contour(X,Y,Z0,[26.8:.1:28.3],'color',[.3 .3 .3]*0);
clabel(c,h,'fontsize',6,'color',[.1 .1 .1])
plot(sa_tmp(:),ct_tmp(:),'.k')
tmps = sa_tmpv(:,indneg(front_ind));tmpt = ct_tmpv(:,indneg(front_ind));
tmpd = dv_sca(:,indneg(front_ind));tmpv = v_tot_tmp(:,indneg(front_ind));
plot(sa_tmpv(:,indneg(front_ind)),ct_tmpv(:,indneg(front_ind)),'-','color',[.5 .5 .5])
scatter(tmps(:),tmpt(:),5,tmpd(:),'filled')
colorbarnew('v',0.01,1,'d_{front} (km)');
colormap(lansey)

xlim([34.9 35.4]);ylim([-1 10])
xlabel('S_A (g kg^{-1})');ylabel('\Theta (\circC)')


axes('position',[2*llx+lx+0.017 .55 lx*0.9 lx*1.1])
hold on;box on
[c,h]=contour(X,Y,Z0,[26.8:.025:28.4],'color',[.8 .8 .8]);
[c,h]=contour(X,Y,Z0,[26.8:.1:28.3],'color',[.3 .3 .3]*0);
clabel(c,h,'fontsize',6,'color',[.1 .1 .1])
plot(sa_tmp(:),ct_tmp(:),'.k')
scatter(tmps(:),tmpt(:),5,tmpv(:),'filled')
colorbarnew('v',0.01,1,'V_g (m s^{-1})');
colormap(lansey)
caxis(v_lim)

xlim([34.9 35.4]);ylim([-1 10])
xlabel('S_A (g kg^{-1})');ylabel('\Theta (\circC)')


eval(['printHR(''fig/TS_MR' num2str(ll) ''')'])





close all;figure('visible',opt_vis)
%%%%%%%%%%%%%%%%%%%% Vcyclo, N2, PV
pv_lim = [-12.25 -9.75];

axes('position',[.5 .08 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;grid on;box on
pcolor(dd_sec,glid_pv.pp(:,1),real(log10(pvA_tmp)));shading interp
[c,h]=contour(dd_sec,glid_pv.pp(:,1),sigs_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(cmocean('curl'));
caxis(pv_lim)
zdn
ylim([0 1000])

ylabel('Depth (m)')
xlabel('Cross-front distance (km)')
colorbarnew('v',0.02,1,'log_{10}(-PV_A)');



axes('position',[.5 0.08+ly+0.05 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;grid on;box on
pcolor(dd_sec,glid_pv.pp(:,1),real(log10(-pvB_tmp)));shading interp
[c,h]=contour(dd_sec,glid_pv.pp(:,1),sigs_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(cmocean('curl'));
caxis(pv_lim)
zdn
ylim([0 1000])

ylabel('Depth (m)')
colorbarnew('v',0.02,1,'log_{10}(-PV_B)');

axes('position',[.5 .08+2*ly+2*0.05 lx ly],'yticklabel',[],'xtick',[0],'xlim',[-100 100],'xticklabel',{datestr(tt_all(ind_front(iq_front)),'dd mmmm, YYYY')},'xaxislocation','top')
axes('position',[.5 0.08+2*ly+2*0.05 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;grid on;box on
pcolor(dd_sec,glid_pv.pp(:,1),real(log10(pv_tmp)));shading interp
[c,h]=contour(dd_sec,glid_pv.pp(:,1),sigs_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(cmocean('curl'));
caxis(pv_lim)
zdn
ylim([0 1000])

ylabel('Depth (m)')
colorbarnew('v',0.02,1,'log_{10}(PV)');

[d_grid p_grid] = meshgrid(dd_sec,glid_pv.pp(:,1));
[phiB phiC iq_stable iq_IISI iq_SI iq_GISI iq_GI] = get_insta_from_PV(M4_tmp,ff_tmp,N2_tmp,vort_tmp);

tit_legend = {'stable','Gravitational','Symmetric/Gravitational','Symmetric','Inertial/Symmetric'};
axes('position',[.07 0.08+2*ly+2*0.05 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;box on;grid on
h1 = plot(d_grid(iq_stable),p_grid(iq_stable),'.r');
h2 = plot(d_grid(iq_GI),p_grid(iq_GI),'.g');
h3 = plot(d_grid(iq_GISI),p_grid(iq_GISI),'.b');
h4 = plot(d_grid(iq_SI),p_grid(iq_SI),'.y');
h5 = plot(d_grid(iq_IISI),p_grid(iq_IISI),'.','color',[.5 .5 .5]);
contour(d_grid,p_grid,[sigs_tmp],[27.3 27.6 27.7 27.74 27.76 27.79:0.02:28.06],'k');
colormap(cmocean('curl')); zdn; ylim([0 1000])
xlabel('Cross-front distance (km)')
ylabel('Depth (m)')

axes('position',[.05 0.08+2*0.05+0.2 lx ly])
hold on;axis off
h1 = plot(d_grid(1),p_grid(1),'.r');h2 = plot(d_grid(1),p_grid(1),'.g');
h3 = plot(d_grid(1),p_grid(1),'.b');h4 = plot(d_grid(1),p_grid(1),'.y');
h5 = plot(d_grid(1),p_grid(1),'.','color',[.5 .5 .5]);
legend([h1 h2 h3 h4 h5],tit_legend,'location','southeast')
axis off



eval(['printHR(''fig/Insta_MR' num2str(ll) ''')'])




%%% watermass distribution
col_paired = paired();
close all;figure('visible',opt_vis)

axes('position',[llx+lx+llx-0.01 llx+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100])
hold on;box on;grid on
h3 = plot(glid_MR.bin_d(glid_MR.maskio_AW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_AW(:,:,ll)),'.','color',col_paired(5,:));
h4 = plot(glid_MR.bin_d(glid_MR.maskio_mAW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_mAW(:,:,ll)),'.','color',col_paired(7,:));
h5 = plot(glid_MR.bin_d(glid_MR.maskio_nAW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_nAW(:,:,ll)),'.','color',col_paired(6,:));
h7 = plot(glid_MR.bin_d(glid_MR.maskio_NSDW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_NSDW(:,:,ll)),'.','color',col_paired(3,:));
h6 = plot(glid_MR.bin_d(glid_MR.maskio_NSIW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_NSIW(:,:,ll)),'.','color',col_paired(4,:));
h1 = plot(glid_MR.bin_d(glid_MR.maskio_sASW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_sASW(:,:,ll)),'.','color',col_paired(1,:));
h2 = plot(glid_MR.bin_d(glid_MR.maskio_wASW(:,:,ll)),glid_MR.bin_p(glid_MR.maskio_wASW(:,:,ll)),'.','color',col_paired(2,:));

[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(cmocean('haline'));zdn;ylim([0 1000])
xlabel('Cross-front distance (km)')

axes('position',[.05 llx-0.05 lx ly])
hold on;axis off
ind_col = [5 7 6 3 4 1 2];tit_legend = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'};
for ii=1:length(ind_col)
eval(['h' num2str(ii) ' = plot(d_grid(1),p_grid(1),''.'',''color'',col_paired(ind_col(ii),:));'])
end
xlim([0 1])
lgd=legend([h1 h2 h3 h4 h5 h6 h7],tit_legend,'location','southeast','NumColumns',7);lgd.Position = [.03 0.02 .9 .03];legend('boxoff')
axis off

ax=axes('position',[llx llx+ly+0.02+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xticklabel','','xlim',[-100 100]);
hold on;box on;grid on
pcolor(dd_sec,0:1010,ct_tmp);shading interp
[c,h]=contour(dd_sec,0:1010,sig_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(ax,cmocean('dense'));caxis(t_lim);zdn;ylim([0 1000])
ylabel('Depth (m)');

a=axes('position',[llx llx+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
[dgrid_tmp pgrid_tmp] = meshgrid(dd_sec,0:5:1001);ii=0;
for watermass = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
ii=ii+1;
eval(['plot(dgrid_tmp(' watermass{1} '_tmp),pgrid_tmp(' watermass{1} '_tmp),''.'',''color'',col_paired(ind_col(ii),:));'])
end
xlabel('Cross-front distance (km)')

[c,h]=contour(dd_sec,0:1010,sig_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(a,cmocean('dense'));zdn;ylim([0 1000])
ylabel('Depth (m)');

axx=axes('position',[llx+lx+llx-0.01 llx+ly+0.02+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xtick',[-200:50:200],'xticklabel','','xlim',[-100 100]);
hold on;box on;grid on
pcolor(glid_MR.bin_d,glid_MR.bin_p,glid_MR.CT_io(:,:,ll));shading interp
[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(axx,cmocean('thermal'));caxis(t_lim);zdn;ylim([0 1000])
colorbarnew('v',0.01,1,'\Theta (\circC)');


axxx=axes('position',[llx llx+2*ly+2*0.02+0.03 lx ly],'yticklabel',[],'xtick',[max(nanmin(dd_sec),-150) min(nanmax(dd_sec),150)],'xlim',[-100 100],'xticklabel',{[datestr(tt_all(ind_front([1 end])),'dd mmmm, YYYY')]},'xaxislocation','top');
axes('position',[llx llx+2*ly+2*0.02+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xticklabel','','xtick',[-200:50:200],'xlim',[-100 100])
hold on;box on;grid on
pcolor(dd_sec,0:1010,sa_tmp);shading interp
[c,h]=contour(dd_sec,0:1010,sig_tmp,[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(axxx,cmocean('haline'));caxis(s_lim);zdn;ylim([0 1000])
ylabel('Depth (m)')

axxxx=axes('position',[llx+lx+llx-0.01 llx+2*ly+2*0.02+0.03 lx ly],'ytick',[0 200 400 600 800 1000],'xticklabel','','xtick',[-200:50:200],'xlim',[-100 100]);
hold on;box on;grid on
pcolor(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SA_io(:,:,ll));shading interp
[c,h]=contour(glid_MR.bin_d,glid_MR.bin_p,glid_MR.SIG_io(:,:,ll),[27.3 27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
colormap(axxxx,cmocean('dense'));caxis(s_lim);zdn;ylim([0 1000])
colorbarnew('v',0.01,1,'S_A (g kg^{-1})');
title('Optimal interpolation')



eval(['printHR(''fig/watermass_MR' num2str(ll) ''')'])