%%%%
d_all = [];t_all = [];ct_all = [];sa_all = [];lon_all = [];lat_all = [];
for ll=[1:length(sec_beg)]
    ind_front = sec_beg(ll):sec_end(ll);
    if rev(ll) == 0
        d_sec = d_along(ind_front)-d_along(ind_front(1));dd_sec = dc_along(ind_front)-dc_along(ind_front(1));
    else
        d_sec = -fliplr(d_along(ind_front)-d_along(ind_front(end)));dd_sec = -fliplr(dc_along(ind_front)-dc_along(ind_front(end)));
        ind_front = fliplr(ind_front);
    end
    ct_tmp = CT_all(:,ind_front);pp_tmp = PP_all(:,ind_front);sa_tmp = SA_all(:,ind_front);
    llg_glid = lg_all(ind_front); llt_glid = lt_all(ind_front); tt_glid = tt_all(ind_front);
    d_sec = d_sec - d_sec(i_front(ll));dd_sec = dd_sec - dd_sec(i_front(ll));

    %%%%% store all variable
    d_all = [dd_sec d_all];t_all = [tt_glid t_all];lon_all = [llg_glid lon_all];lat_all = [llt_glid lat_all];
    ct_all = [ct_tmp(1:5:end,:) ct_all];sa_all = [sa_tmp(1:5:end,:) sa_all];
end

%%%%% intro map
a = col_paired;
lon_g = [-5.5 5.5];lat_g = [70+2/3 73+1/3];
lg_l = [-15 20];lt_l = [66 80];
glid_l = mergestruct(sg560,sg562);

plot_circulation_LB

close all

ax1 =axes('position',[-.04 .52 .5 .35],'xaxislocation','bottom','yaxislocation','left');
hold on;grid on;box on

m_proj('lambert','lat',lt_l,'lon',lg_l);

m_grid('xtick',[-15:10:20],'tickdir','out','ytick',[65:5:80],'linest','-');

aa=m_pcolor(lg_bat,lt_bat,z);shading interp;set(aa,'facealpha',.7);
m_contour(lg_bat,lt_bat,z_smo,[-4000:500:-500],'color',[.5 .5 .5],'linewidth',.5);

m_gshhs_i('patch',[.6 .7 .6]);
m_gshhs_i('speckle','color','k');
m_text(15,67,'\bfNorway','rotation',0,'fontsize',6,'color','k')

m_plot(lg_deep(1:end),lt_deep(1:end),':','linewidth',2,'color',a(6,:))
m_plot(lg_deep(1:end),lt_deep(1:end),':','linewidth',.5,'color','w')
plot_m_triangle(lg_deep(1:end),lt_deep(1:end),1,a(6,:))

m_plot(lg_west,lt_west,'-','linewidth',2,'color',a(6,:))
m_plot(lg_west,lt_west,'-','linewidth',.5,'color','w')
m_plot(lg_MR(1:end),lt_MR(1:end),'-','linewidth',2,'color',a(6,:))
m_plot(lg_MR(1:end),lt_MR(1:end),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_MR(1:end),lt_MR(1:end),1.5,a(6,:))

m_plot(lg_outfl(1:end-1),lt_outfl(1:end-1),'-','linewidth',2,'color',a(6,:))
m_plot(lg_outfl(1:end-1),lt_outfl(1:end-1),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_outfl(1:end-1),lt_outfl(1:end-1),2,a(6,:))
m_plot(lg_outfl2(1:end),lt_outfl2(1:end),'-','linewidth',2,'color',a(6,:))
m_plot(lg_outfl2(1:end),lt_outfl2(1:end),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_outfl2(1:end),lt_outfl2(1:end),1.5,a(6,:))
m_plot(lg_outfl3(1:end),lt_outfl3(1:end),'-','linewidth',2,'color',a(6,:))
m_plot(lg_outfl3(1:end),lt_outfl3(1:end),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_outfl3(1:end),lt_outfl3(1:end),1.5,a(6,:))
m_plot(lg_recirc(1:end),lt_recirc(1:end),':','linewidth',2,'color',a(5,:))
m_plot(lg_recirc(1:end),lt_recirc(1:end),':','linewidth',.5,'color','w')
plot_m_triangle(lg_recirc(1:end),lt_recirc(1:end),1.5,a(5,:))


m_plot(lg_slope(1:end),lt_slope(1:end),'-','linewidth',2,'color',a(6,:))
m_plot(lg_bso2(1:end-10),lt_bso2(1:end-10),'-','linewidth',2,'color',a(6,:))
m_plot(lg_bso2(1:end-10),lt_bso2(1:end-10),'-','linewidth',.5,'color','w')
m_plot(lg_slope(1:end),lt_slope(1:end),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_bso2(1:end-15),lt_bso2(1:end-15),1,a(6,:))
plot_m_triangle(lg_slope(1:end),lt_slope(1:end),1.5,a(6,:))

m_plot(lg_JMC(25:end),lt_JMC(25:end),'-','linewidth',2,'color',a(1,:))
m_plot(lg_JMC(25:end),lt_JMC(25:end),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_JMC(1:end-5),lt_JMC(1:end-5),1,a(1,:))


m_plot(lg_egc(1:end-20),lt_egc(1:end-20),'-','linewidth',2,'color',a(2,:))
m_plot(lg_egc(1:end-20),lt_egc(1:end-20),'-','linewidth',.5,'color','w')
plot_m_triangle(lg_egc(1:end-20),lt_egc(1:end-20),1.5,a(2,:))

m_plot([lon_g fliplr(lon_g) lon_g(1)],[lat_g(1) lat_g fliplr(lat_g)],'-k','linewidth',.5)

m_text(-13,78,'\bfEGC','rotation',0,'fontsize',7,'color',a(2,:))
m_text(7,66.75,'\bfNwASC','rotation',0,'fontsize',7,'color',a(6,:))
m_text(-3,66.5,'\bfNwAFC','rotation',0,'fontsize',7,'color',a(6,:))
m_text(12,78,'\bfWSC','rotation',0,'fontsize',7,'color',a(6,:))
m_text(-4,75,'\bfGS','rotation',0,'fontsize',7,'color','k')
m_text(2,70,'\bfLB','rotation',0,'fontsize',7,'color','k')
m_text(4,67,'\bfVP','rotation',0,'fontsize',7,'color','k')
m_text(0,72,'\bfMR','rotation',0,'fontsize',7,'color','k')
m_text(7,74,'\bfKR','rotation',0,'fontsize',7,'color','k')
m_text(-21,72,'\bfJMC','rotation',0,'fontsize',7,'color',a(1,:))

caxis([-3500 -1000]);h=colorbarnew('h',-0.12,0.5,'Bathymetry (m)');
cbpos = get(h,'position');set(h,'position',[cbpos(1)+.055 cbpos(2)+.23 cbpos(3) cbpos(4)]);
set(h,'xtick',[-3500:1000:-1500],'xaxislocation','top');

labfig('\bfa',-2.8,-0.3,0);

%% subplot 2
m_proj('equidist','lat',lat_g,'lon',lon_g);

ax2 = axes('position',[.38 .43 .52 .57],'yaxislocation','right');
hold on;grid on;box on
m_grid('linestyle','none','tickdir','out','ytick',66:1/3:74);

aa=m_pcolor(lg_bat,lt_bat,z);shading interp;set(aa,'facealpha',.7);
quiv_gl = m_quiver([lon_g(1)+0.3 [glid_l.LGV]],[lat_g(2)-0.52 [glid_l.LTV]],[0.3 [glid_l.VX]],[0 [glid_l.VY]],2,'color',[.2 .2 .2],'showarrowhead','off');
for l=1:2
plt_glid(l) = m_plot(nanmean(glid_l(l).ZLG),nanmean(glid_l(l).ZLT),'-','color',[.2 .2 .2]);
end
m_text(lon_g(1)+0.3,lat_g(2)-0.45,'0.3 m s^{-1}','fontsize',6)

%%mean section
m_plot(lggg,lttt,'-k','linewidth',1.5);m_plot(lggg,lttt,'-w','linewidth',0.5);
m_plot(nanmean(lg_front),nanmean(lt_front),'ok','markerfacecolor','k','markersize',3.5)
m_plot(nanmean(lg_front),nanmean(lt_front),'ow','markerfacecolor','w','markersize',2.5)
%%%% reference frame
aa=(lttt(2)-lttt(1))/(lggg(2)-lggg(1));lgO = -1.95;ltO = 73.12;
m_plot(lgO+[0 -.25],ltO+[0 -aa/4],'-k','linewidth',1.5);m_plot(lgO-[0 aa/4]/cosd(ltO),ltO-[0 -.25]*cosd(ltO),'-k','linewidth',1.5)
plot_m_triangle(lgO+[0 -.25],ltO+[0 -aa/4],.15,[0 0 0]);plot_m_triangle(lgO-[0 aa/4]/cosd(ltO),ltO-[0 -.25]*cosd(ltO),.15,[0 0 0])
m_text(lgO-.25-.3,ltO-aa/4,'y','fontsize',8);m_text(lgO-aa/4/cosd(ltO)+.1,ltO+.25*cosd(ltO)-.05,'x','fontsize',8);

% transport
h1=m_quiver([lon_g(1)+0.3 lg_front],[lat_g(2)-0.37 lt_front],[5 dir_front(1,:).*TrAW],[0 dir_front(2,:).*TrAW],.5,'color','w','ShowArrowHead','off','linewidth',2);
h1=m_quiver([lon_g(1)+0.3 lg_front],[lat_g(2)-0.37 lt_front],[5 dir_front(1,:).*TrAW],[0 dir_front(2,:).*TrAW],.5,'color','r','ShowArrowHead','off','linewidth',1);

m_plot(lg_front,lt_front,'ok','markerfacecolor','w','markersize',7)

for ll=1:length(lg_front)
if ll<10
m_text(lg_front(ll)-.07,lt_front(ll),['\bf' num2str(ll)],'fontsize',5)
elseif ll==10
m_text(lg_front(ll)-.07,lt_front(ll)+.2,['\bf' num2str(ll)],'fontsize',5)
m_plot(lg_front(ll)+[0 0.07],lt_front(ll)+[0.06 .16],'-k')
else
m_text(lg_front(ll)-.13,lt_front(ll),['\bf' num2str(ll)],'fontsize',5)
end
end
m_text(lon_g(1)+0.3,lat_g(2)-0.3,'5 Sv','fontsize',6)

caxis([-3500 -1000]);
h=colorbarnew('h',0,0.3,'Depth (m)');

cbpos = get(h,'position');set(h,'position',[cbpos(1)+.34 cbpos(2)-0.01 cbpos(3) cbpos(4)],'color','k');
m_ruler([.02 .22]+.76,.95,3,'ticklen',.01,'fontsize',5);

labfig('\bfb',-0.2,-0.4,0);



%%%% hovmueller plot!
ax3 = axes('position',[.1 .21 .8 .14],'yaxislocation','left');
hold on;grid on;box on

ind = find(bin_p>0 & bin_p<50);
m_plot(tt_front,zeros(size(tt_front)),'ok','markerfacecolor','k','markersize',4)

plot(t_all,d_all,'ok','markersize',3)
scatter(t_all,d_all,6,nanmean(sa_all(ind,:)),'filled')

for ll=1:length(tt_front)
text(tt_front(ll)-2-2*(ll>=10),165,num2str(ll),'fontsize',7)
end
text(datenum(2016,6,5),-90,'\bf0-25m','fontsize',6)

ylim([-110 110]);caxis([34.95 35.35])
xlim([datenum(2016,5,27) datenum(2017,8,1)]);set(gca,'xtick',[datenum(2016,5:12,1) datenum(2017,1:8,1)])
nolabx

labfig('\bfc',-0.05,-1,0);

ylabel('y (km)')
h=colorbarnew('v',0.005,1,'S_A (g kg^{-1})');


ax4 = axes('position',[.1 .06 .8 .14],'yaxislocation','left');
hold on;grid on;box on

% project onto mean section
[xx_all yy_all] = earthcoord('flat',lg_all,lt_all,nanmean(lg_front),nanmean(lt_front));
[x_p y_p] = project_pt(xx_all,yy_all,-nanmean(dir_front(1,:))/nanmean(dir_front(2,:)),0);

ind = find(bin_p>300 & bin_p<400);

m_plot(tt_front,zeros(size(tt_front)),'ok','markerfacecolor','k','markersize',4)

plot(t_all,d_all,'ok','markersize',3)
scatter(t_all,d_all,6,nanmean(ct_all(ind,:)),'filled')
text(datenum(2016,6,5),-90,'\bf300-400m','fontsize',6)

ylim([-110 110]);caxis([-0.5 5.5])
xlim([datenum(2016,5,27) datenum(2017,8,1)]);set(gca,'xtick',[datenum(2016,5:12,1) datenum(2017,1:8,1)])
dtickx


labfig('\bfd',-0.05,-1,0);

ylabel('y (km)')
h=colorbarnew('v',0.005,1,'\Theta (\circC)');



axes('position',[.1 .06-0.05 .8 .02])
hold on; axis off
text(datenum(2016,9,1),0,'2016','fontsize',6)
text(datenum(2017,1,1),0,'|','fontsize',6)
text(datenum(2017,4,1),0,'2017','fontsize',6)
xlim([datenum(2016,5,25) datenum(2017,8,1)])

colormap(ax1,cmocean('-deep'));colormap(ax2,cmocean('-deep'));colormap(ax3,cmocean('haline'));colormap(ax4,cmocean('thermal'));

printHR('fig/map_MR_intro')