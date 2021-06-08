%%%%%% not sure beta effect is imporant here...

[glid_MR.maskio_sASW glid_MR.maskio_wASW glid_MR.maskio_AW glid_MR.maskio_mAW glid_MR.maskio_nAW glid_MR.maskio_NSIW glid_MR.maskio_NSDW] = select_water_mass(glid_MR.SA_io,glid_MR.CT_io);
glid_MR.maskio_ASW = ceil((glid_MR.maskio_wASW+glid_MR.maskio_sASW)/2);

close all
lx = .24; ly = .25;

ax1=axes('position',[.08 .7 lx ly]);
hold on;grid on;box on
aa=pcolor(dd_grid,pp_grid,ct_mean);shading interp;set(aa,'facealpha',.9);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.5);
clabel(c,h,'fontsize',5,'color',[.3 .3 .3]);
zdn;xlim([-85 85]);
caxis([-1.5 8])
h=colorbarnew('v',0.005,1);set(get(h,'title'),'string','\Theta (\circC)')
ylabel('Depth (m)')
set(gca,'xlim',[-85 85],'ylim',[0 1000],'ytick',0:250:1000,'xtick',-90:30:90)
labfig('\bfa',-0.25,-0.75,1);
xlabel('y (km)')

ax2=axes('position',[.08+lx+0.07 .7 lx ly]);
hold on;grid on;box on
aa=pcolor(dd_grid,pp_grid,sa_mean);shading interp;set(aa,'facealpha',.9);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.5);
clabel(c,h,'fontsize',5,'color',[.3 .3 .3]);
zdn;xlim([-85 85]);nolaby
caxis([35.05 35.35])
h=colorbarnew('v',0.005,1);set(get(h,'title'),'string','S_A (g kg^{-1})')
set(gca,'xlim',[-85 85],'ylim',[0 1000],'ytick',0:250:1000,'xtick',-90:30:90)
labfig('\bfb',-0.25,-0.75,1);
xlabel('y (km)')


ax3=axes('position',[.08+2*(lx+0.07) .7 lx ly]);
hold on;grid on;box on
aa=pcolor(dd_grid,pp_grid,-vg_mean);shading interp;set(aa,'facealpha',.9);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(dd_grid,pp_grid,-vg_mean,-[-0.55:.1:.15],'-','color',[.3 .3 .3]);
clabel(c,h,'fontsize',5,'color',[.3 .3 .3]);
zdn;xlim([-85 85]);nolaby
caxis([-.45 .45])
h=colorbarnew('v',0.005,1);set(get(h,'title'),'string','U_g (m s^{-1})');
set(gca,'xlim',[-85 85],'ylim',[0 1000],'ytick',0:250:1000,'xtick',-90:30:90)
labfig('\bfc',-0.25,-0.75,1);
xlabel('y (km)')


colormap(ax1,cmocean('thermal'));colormap(ax2,cmocean('haline'));colormap(ax3,cmocean('curl'));

printHR('fig/paper_TSV')


%%%% seasons
lab_l = {'\bfa','\bfb','\bfc','\bfd','\bfe','\bff','\bfg','\bfh','\bfi','\bfj','\bfk','\bfl'};
seas_l = {'JFMA','MJJA','SOND'};
tit_l = {'\Theta (\circC)','S_A (g kg^{-1})','U_g (m s^{-1})','Obs. Probability'};
ax_l = {[-0.5 7],[35.05 35.35],[-.45 .45],[0 3]};
ncol = 4;nlig = 3;
close all
[ax,row,col,IsTop,IsBottom,IsLeft,IsRight]=subaxes(nlig,ncol,'RowSpacing',0.1,'ColSpacing',0.1,'PageMargins',[0.08 0.09 0.08 0.1]);
for l=1:length(ax)
hold(ax(l),'on');grid(ax(l),'on');box(ax(l),'on');set(ax(l),'ydir','rev');
set(ax(l),'xlim',[-85 85],'ylim',[0 1000],'xtick',-90:30:90)
if ~IsLeft(l)
set(ax(l),'yticklabel','')
else
ylabel(ax(l),'Depth (m)')
end
if ~IsBottom(l)
set(ax(l),'xticklabel','')
else
xlabel(ax(l),'y (km)')
end
if IsTop(l)
eval(['cb' num2str(col(l)) '=colorbarnew(''h'',-.25,1,tit_l{l},ax(l));']) % set(get(cb,'title'),'String',tit_l{l})
if col(l) == 3
eval(['set(cb' num2str(col(l)) ',''ylim'',[-.45 .15]);'])
end
end
end


%%%% seasonal
varr = {'ct','sa','vg'};
iq = find((row == 1 | row == 2 | row == 3) & (col == 1 | col == 2 | col == 3));
for ll=iq
eval(['vartmp = ' varr{col(ll)} '_seas(:,:,row(ll));'])
aa=pcolor(ax(ll),dd_grid,pp_grid,vartmp);shading(ax(ll),'interp');set(aa,'facealpha',.9);
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_seas(:,:,row(ll)),[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
if col(ll)==3
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_seas(:,:,row(ll)),[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(ax(ll),dd_grid,pp_grid,vartmp,[-.45:0.1:.15],'color',[.3 .3 .3],'linewidth',.5);
clabel(c,h,'fontsize',4,'color',[.3 .3 .3]);
else
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_seas(:,:,row(ll)),[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.5);
end
if col(ll)==1
clabel(c,h,'fontsize',4,'color',[.3 .3 .3]);
a=text(60,900,['\bf' seas_l{row(ll)}],'fontsize',5);set(a,'parent',ax(ll));
end
caxis(ax(ll),ax_l{col(ll)})
labfig(lab_l{ll},-0.25,-0.75,1,ax(ll));
end


iq = find(row == 4 & col<4);
for ll=iq
eval(['vartmp = ' varr{col(ll)} '_mean;'])
pcolor(ax(ll),dd_grid,pp_grid,vartmp);shading(ax(ll),'interp');
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_mean,[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_mean,[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.1);
if col(ll)==1
clabel(c,h,'fontsize',4,'color',[.3 .3 .3]);
a=text(60,900,'\bfannual','fontsize',5);set(a,'parent',ax(ll));
end
caxis(ax(ll),ax_l{col(ll)})
end

iq = find(col == 1);
for ll=iq
colormap(ax(ll),cmocean('thermal'))
end
iq = find(col == 2);
for ll=iq
colormap(ax(ll),cmocean('haline'))
end
iq = find(col == 3);
for ll=iq
colormap(ax(ll),cmocean('curl'))
end
%%%% annula av
load cmap_proba
col_samp = [14 5 1];cmap_proba = [];
for ii=1:length(col_samp)
cmap_proba = [cmap_proba ; cmap_pb(col_samp(ii)*16+[1:8],:)];
end

seas = {'JFMA','MJJA','SOND'};iq = find(col == 4);
for ll=iq
eval(['iq_seas = iq' seas{row(ll)} ';'])
for watermass = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
eval(['p' watermass{1} '_seas(:,:,row(ll)) = nansum(glid_MR.maskio_' watermass{1} '(:,:,iq_seas),3)./nobs_seas(:,:,row(ll));p' watermass{1} '_seas(p' watermass{1} '_seas==0|p' watermass{1} '_seas>1|nobs_seas(:,:,row(ll))<4)=NaN;'])
end
pASW_seas(:,:,row(ll))=nansum(glid_MR.maskio_sASW(:,:,iq_seas)+glid_MR.maskio_wASW(:,:,iq_seas),3)./nobs_seas(:,:,row(ll));pASW_seas(pASW_seas==0|pASW_seas>1)=NaN;
pDW_seas(:,:,row(ll))=nansum(glid_MR.maskio_NSIW(:,:,iq_seas)+glid_MR.maskio_NSDW(:,:,iq_seas),3)./nobs_seas(:,:,row(ll));pDW_seas(pDW_seas==0|pDW_seas>1)=NaN;
aa=pcolor(ax(ll),dd_grid,pp_grid,1+pAW_seas(:,:,row(ll))*.99);shading(ax(ll),'interp');set(aa,'facealpha',.8);
aa=pcolor(ax(ll),dd_grid,pp_grid,2+pDW_seas(:,:,row(ll))*.99);shading(ax(ll),'interp');set(aa,'facealpha',.8);
aa=pcolor(ax(ll),dd_grid,pp_grid,pASW_seas(:,:,row(ll))*.99);shading(ax(ll),'interp');set(aa,'facealpha',.8);
[c,h]=contour(ax(ll),dd_grid,pp_grid,vg_seas(:,:,row(ll)),[-.05 -.05],'--','color',[0 0 0],'linewidth',.75);
contour(ax(ll),dd_grid,pp_grid,~isnan(pmAW_seas(:,:,row(ll))),[0.5 0.5],'color',col_paired(7,:),'linewidth',.75);
%contour(ax(ll),dd_grid,pp_grid,~isnan(pmAW_seas(:,:,row(ll))),[0.5 0.5],'color','w','linewidth',.1);
contour(ax(ll),dd_grid,pp_grid,~isnan(pnAW_seas(:,:,row(ll))),[0.5 0.5],'color',col_paired(5,:),'linewidth',.75);
%contour(ax(ll),dd_grid,pp_grid,~isnan(pnAW_seas(:,:,row(ll))),[0.5 0.5],'color','w','linewidth',.1);
%contour(ax(ll),dd_grid,pp_grid,~isnan(pDW_seas(:,:,row(ll))),[0.5 0.5],'color',col_paired(10,:));
contour(ax(ll),dd_grid,pp_grid,~isnan(pASW_seas(:,:,row(ll))),[0.5 0.5],'color',col_paired(1,:),'linewidth',.75);
%contour(ax(ll),dd_grid,pp_grid,~isnan(pASW_seas(:,:,row(ll))),[0.5 0.5],'color','w','linewidth',.1);
colormap(ax(ll),cmap_proba);
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_seas(:,:,row(ll)),[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(ax(ll),dd_grid,pp_grid,sig_seas(:,:,row(ll)),[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.1);
caxis(ax(ll),ax_l{col(ll)})
if row(ll)==1
aa=text(-85,450,'mAW','color',col_paired(7,:));set(aa,'parent',ax(ll),'fontsize',7)
elseif row(ll)==2
aa=text(-50,75,'nAW','color',col_paired(5,:));set(aa,'parent',ax(ll),'fontsize',7)
end
labfig(lab_l{ll},-0.25,-0.75,1,ax(ll));
end
set(cb4,'ytick',[0:.5:3],'yticklabel',{'0','ASW','1','AW','1','DW','1'})

printHR('fig/paper_TSV_seas')





%%% new figure about water mass
lx = .3; ly = .45;

close all
ax7=axes('position',[.08 .5 lx ly]);
hold on;grid on;box on

iqMR = find(lg_all<4.5 & lt_all>70.75);
glidt = CT_all(:,iqMR);glids = SA_all(:,iqMR);

[mz,xmid,ymid,numz,stdz]=twodstats([glids(:)'],[glidt(:)'],ones(size([glids(:)'])),200);
numz(numz==0)=NaN;
[sgrid tgrid] = meshgrid(xmid,ymid);
sig0 = gsw_sigma0(sgrid,tgrid);sig1 = gsw_sigma1(sgrid,tgrid);

%% water mass
aa=pcolor(sgrid,tgrid,log10(numz));shading flat;set(aa,'facealpha',.9)
[c,h] = contour(sgrid,tgrid,sig0,[26:.1:28.15],'color',[.3 .3 .3]);
clabel(c,h,'fontsize',4,'color',[.3 .3 .3]);
caxis([-1 4])
h=colorbarnew('v',0.005,1);set(get(h,'title'),'string','# count');set(h,'ylim',[0 4]);set(h,'ytick',[0:4],'yticklabel',{'1';'10';'10^2';'10^3';'10^4'});
xlabel('S_A (g kg^{-1})');ylabel('\Theta (\circC)')
%% labels
[i_sASW i_wASW i_AW i_mAW i_nAW i_NSIW i_NSDW] = select_water_mass(sgrid,tgrid);
ind_col = [6 7 5 7 8 1 2];tit_legend = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'};
for ii=1:length(ind_col)
eval(['contour(sgrid,tgrid,i_' tit_legend{ii} ',''-'',''color'',col_paired(ind_col(ii),:),''linewidth'',0.5);'])
end
text(35.18,10,'\bfAW','color',col_paired(6,:),'fontsize',7);
text(35.33,10,'nAW','color',col_paired(5,:),'fontsize',7);
text(35.3,4.5,'mAW','color',col_paired(7,:),'fontsize',7);
text(35.1,-0.5,'NSDW','color',col_paired(7,:),'fontsize',7);
text(34.97,-0.5,'NSIW','color',col_paired(8,:),'fontsize',7);
text(34.93,1.5,'wASW','color',col_paired(2,:),'fontsize',7);
text(35.05,7.5,'sASW','color',col_paired(1,:),'fontsize',7);
labfig('\bfa',1,1.5,0);

load stern.mat
colormap(ax7,stern);

%ax3=axes('position',[.1 .1 .8 .3])
%plot(month_front,Tr_AW,'o')

for watermass = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
eval(['p' watermass{1} ' = nansum(glid_MR.maskio_' watermass{1} ',3)./nobs;p' watermass{1} '(p' watermass{1} '==0|p' watermass{1} '>1)=NaN;'])
end
pASW=nansum(glid_MR.maskio_sASW+glid_MR.maskio_wASW,3)./nobs;pASW(pASW==0|pASW>1)=NaN;
pDW=nansum(glid_MR.maskio_NSIW+glid_MR.maskio_NSDW,3)./nobs;pDW(pDW==0|pDW>1)=NaN;
ax8=axes('position',[.08+lx+0.15 .5 .4 ly]);
hold on;grid on;box on
aa=pcolor(dd_grid,pp_grid,pAW*.99);shading flat;set(aa,'facealpha',.7);
%aa=pcolor(dd_grid,pp_grid,1+pnAW);shading interp;set(aa,'facealpha',.5)
%aa=pcolor(dd_grid,pp_grid,2+pmAW);shading interp;set(aa,'facealpha',.5)
aa=pcolor(dd_grid,pp_grid,1+pDW*.99);shading flat;set(aa,'facealpha',.7);
aa=pcolor(dd_grid,pp_grid,2+pASW*.99);shading flat;set(aa,'facealpha',.9);
contour(dd_grid,pp_grid,~isnan(pAW),[0.5 0.5],'color',col_paired(6,:));
contour(dd_grid,pp_grid,~isnan(pDW),[0.5 0.5],'color',col_paired(8,:));
contour(dd_grid,pp_grid,~isnan(pASW),[0.5 0.5],'color',col_paired(2,:));
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.2:.1:27.6 27.8:0.05:28.15],'color',[.3 .3 .3],'linewidth',.1);
[c,h]=contour(dd_grid,pp_grid,sig_mean,[27.3:.2:27.7 27.75:0.1:28.1],'color',[.3 .3 .3],'linewidth',.1);
zdn;xlim([-85 85])
h=colorbarnew('v',0.005,1);set(get(h,'title'),'string','Probability')
ylabel('Depth (m)')
set(gca,'xlim',[-85 85],'ylim',[0 1000],'ytick',0:250:1000,'xtick',-90:30:90)
xlabel('y (km)')
labfig('\bfb',1,1.5,1);
colormap(ax8,cmap_proba);caxis([0 3]);set(h,'yticklabel',[0 0.5 1 0.5 1 0.5 1])

printHR('fig/paper_water_seas')