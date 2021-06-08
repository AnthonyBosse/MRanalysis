close all
clear all
warning off
addpath(genpath('.'))
def_fontname % my customized font for plots

% check if data are already in data folder
clear sg560 sg562
z = dir('data');
if isempty(find(strcmp({z.name},'final_sg560_MR1.nc'))) | isempty(find(strcmp({z.name},'final_sg562_MR2.nc')))
    ftpobj = ftp('ftp.nmdc.no');
    cd(ftpobj,'nmdc/UIB/PROVOLO/Glider/');
    dir(ftpobj)
    if isempty(find(strcmp({z.name},'final_sg560_MR1.nc'))); mget(ftpobj,'final_sg560_MR1.nc','data'); end
    if isempty(find(strcmp({z.name},'final_sg562_MR2.nc'))); mget(ftpobj,'final_sg562_MR2.nc','data'); end
end
% load data
sg560.glider_name = 'sg560';sg562.glider_name = 'sg562';
sg560.deployment = 'MR1';sg562.deployment = 'MR2';
for file = {'final_sg560_MR1','final_sg562_MR2'}
        l=0;var_short={'JULD','LG','LT','PP','T','SS','DAC_JULD','LGDAC','LTDAC','U','V'};
        for varr = {'JULD','LONGITUDE','LATITUDE','PRES','TEMP','PSAL','DAC_JULD','DAC_LONGITUDE','DAC_LATITUDE','DAC_U','DAC_V'}
        l = l+1;
        eval([file{1}(7:11) '.' var_short{l} ' = ncread(''' file{1} '.nc'',''' varr{1} ''');'])
    end
end
sg560.DAYS = sg560.JULD + datenum(1950,01,01,00,00,00);sg560.DAC_TIME = sg560.DAC_JULD + datenum(1950,01,01,00,00,00);
sg562.DAYS = sg562.JULD + datenum(1950,01,01,00,00,00);sg562.DAC_TIME = sg562.DAC_JULD + datenum(1950,01,01,00,00,00);

% interpolate DAC estimate on profiles
sg560.VX = interp1(sg560.DAC_TIME,sg560.U,sg560.DAYS)';sg562.VX = interp1(sg562.DAC_TIME,sg562.U,sg562.DAYS)';
sg560.VY = interp1(sg560.DAC_TIME,sg560.V,sg560.DAYS)';sg562.VY = interp1(sg562.DAC_TIME,sg562.V,sg562.DAYS)';
sg560.TVEL = sg560.DAYS';sg562.TVEL = sg562.DAYS';sg560.LGV = sg560.LG';sg562.LGV = sg562.LG';sg560.LTV = sg560.LT';sg562.LTV = sg562.LT';
sg560.ZDIVE = nansum(~isnan(sg560.T));sg562.ZDIVE = nansum(~isnan(sg562.T));

% get Conservative Temperature, Absolute Salinity and Potential Density matrix
sg560.SA = gsw_SA_from_SP(sg560.SS,sg560.PP,sg560.LG,sg560.LT);sg562.SA = gsw_SA_from_SP(sg562.SS,sg562.PP,sg562.LG,sg562.LT);
sg560.CT = gsw_CT_from_t(sg560.SA,sg560.T,sg560.PP);sg562.CT = gsw_CT_from_t(sg562.SA,sg562.T,sg562.PP);
sg560.SIG = gsw_sigma0(sg560.SA,sg560.CT);sg562.SIG = gsw_sigma0(sg562.SA,sg562.CT);

% expand time, longitude, latitude matrices in 2D (assuming vertical profiles here)
sg560.DAYS = repmat(sg560.DAYS',size(sg560.PP,1),1);sg562.DAYS = repmat(sg562.DAYS',size(sg562.PP,1),1);
sg560.ZLG = repmat(sg560.LG',size(sg560.PP,1),1);sg562.ZLG = repmat(sg562.LG',size(sg562.PP,1),1);
sg560.ZLT = repmat(sg560.LT',size(sg560.PP,1),1);sg562.ZLT = repmat(sg562.LT',size(sg562.PP,1),1);
sg560.PP = repmat(sg560.PP,1,size(sg560.DAYS,2));sg562.PP = repmat(sg562.PP,1,size(sg562.DAYS,2));

sg560.VX(sg560.ZDIVE<700) = NaN;sg560.VY(sg560.ZDIVE<700) = NaN;
sg562.VX(sg562.ZDIVE<700) = NaN;sg562.VY(sg562.ZDIVE<700) = NaN;

sg562

% compute Vg and other PV diagnosticis (this step might take a bit of time)
[glid_pv1] = comp_runningPV(sg560,5);[glid_pv2] = comp_runningPV(sg562,5);

% Merge glid_pv1 and glid_pv2
varr = fieldnames(glid_pv1);
for ll=3:length(varr)
    eval(['glid_pv.' varr{ll} ' = [glid_pv1.' varr{ll} ' glid_pv2.' varr{ll} '];'])
end
glid_pv.d_along = [glid_pv1.d_along  glid_pv1.d_along(end)+glid_pv2.d_along];
glid_pv.d_cor = nancumsum([0 diff(nanmean(glid_pv.d_along)).*red_x(nanmean(glid_pv.sint))]);

glid_pv

% some colorbar limits
t_lim = [-1 8];s_lim = [35.05 35.4];sp_lim = [0 1];ox_lim = [295 345];
lon_MR = [-3 6];lat_MR = [70.5 73];
wpt_lt = [71+42/60 72+48/60 72+27.69/60 71+22.56/60 71+1.38/60 72+7.32/60 71+47.65/60 70+49.54/60 71+7.44/60 72];
wpt_lg  = [3+18/60 42/60 0.66/60 2+25.05/60 1+22.596/60 -1-4.488/60 -2-7.08/60 0 14.73/60 2];

% load bathy
load Nordic_ETOPO2.mat
lg_bat = lon;lt_bat = lat;z = elev;
lg_l = [-17 23];lt_l = [61 80];
dx = 25;
indlg_in = find(lg_bat>lg_l(1) & lg_bat<lg_l(2));lg_bat = lg_bat(indlg_in);
indlt_in = find(lt_bat>lt_l(1) & lt_bat<lt_l(2));lt_bat = lt_bat(indlt_in);
z = z(indlt_in,indlg_in);
[lgg_bat ltt_bat] = meshgrid(lg_bat,lt_bat);

% merge DACs and T,S,SIG
TVEL_all = [sg560.TVEL sg562.TVEL];VX_all = [sg560.VX sg562.VX];VY_all = [sg560.VY sg562.VY];
LGV_all = [sg560.LGV sg562.LGV];LTV_all = [sg560.LTV sg562.LTV];
lg_all = nanmean([sg560.ZLG sg562.ZLG]);lt_all = nanmean([sg560.ZLT sg562.ZLT]);tt_all = nanmean([sg560.DAYS sg562.DAYS]);
SA_all = [sg560.SA sg562.SA];CT_all = [sg560.CT sg562.CT];SIG_all = [sg560.SIG sg562.SIG];PP_all = [sg560.PP sg562.PP];DAYS_all = [sg560.DAYS sg562.DAYS];

% compute MLD with 0.03 kg m^-3 criterion
for l=1:size(CT_all,2)
    MLD_all(l) = compute_mld(PP_all(:,l),SIG_all(:,l),0.03);
end
% get surface, Atlantic Water, deep characteristics
ind_surf = find(nanmean(SA_all(1:20,:))-35.05<0); tmpp = diff(ind_surf); iq_surffront = [ind_surf(1) ind_surf(find(tmpp>3)+1)];
ind_atl = find(nanmean(CT_all(300:400,:))-4>0); tmpp = diff(ind_atl); iq_atlfront = [ind_atl(1) ind_atl(find(tmpp>1)+1)];
ind_deep = find(nanmean(SIG_all(950:end,:))>28.06); tmpp = diff(ind_deep); iq_deepfront = [ind_deep(1) ind_deep(find(tmpp>3)+1)];

% divide the stdudy period in 4 periods
ind1 = 40:465; ind2 = 466:1171; ind3 = 1407:1930;ind4 = 1930:2440;
ind1v = find(TVEL_all>tt_all(ind1(1)) & TVEL_all<tt_all(ind1(end)));
ind2v = find(TVEL_all>tt_all(ind2(1)) & TVEL_all<tt_all(ind2(end)));
ind3v = find(TVEL_all>tt_all(ind3(1)) & TVEL_all<tt_all(ind3(end)));
ind4v = find(TVEL_all>tt_all(ind4(1)) & TVEL_all<tt_all(ind3(end)));

close all
lx = .195;ly =.4;
m_proj('equidist','lat',lat_MR,'lon',lon_MR);
for l=1:4
eval(['ind = ind' num2str(l) ';'])
%%%%%%

ax1=axes('position',[.05+(lx+0.02)*(l-1) .5 lx ly]);
if l~=1
    m_grid('linestyle','none','tickdir','out','yticklabel',[])
else
    m_grid('linestyle','none','tickdir','out')
end
hold on;grid on;box on
m_contour(lg_bat,lt_bat,z,[-4000:500:0],'color',[.5 .5 .5])
m_quiver([lon_MR(2)-2 lg_all(ind)],[lat_MR(1)+0.1 lt_all(ind)],[.2 abs(nanmean(glid_pv.v_tot_cor(1:100,ind))).*glid_pv.dir_flow(1,ind)],[0 abs(nanmean(glid_pv.v_tot_cor(1:100,ind))).*glid_pv.dir_flow(2,ind)]*cosd(nanmean(lat_MR)),5,'color',[.2 .2 .2]*0,'ShowArrowHead','off')

m_scatter(lg_all(ind),lt_all(ind),5,nanmean(CT_all(300:400,ind)),'filled')
m_text(lon_MR(2)-2,lat_MR(1)+0.3,'0.2 m s^{-1}','fontsize',5)
colormap(ax1,cmocean('thermal'))

caxis([-0.5 5.5])
if l==4
    h=colorbarnew('v',0.02,1,'\Theta^{300-400} (\circC)');
end
title(['\bf' datestr(tt_all(ind(1)),'mmm YYYY') ' - ' datestr(tt_all(ind(end)),'mmm YYYY')])


ax2=axes('position',[.05+(lx+0.02)*(l-1) .07 lx ly]);
if l~=1
    m_grid('linestyle','none','tickdir','out','yticklabel',[])
else
    m_grid('linestyle','none','tickdir','out')
end
hold on;grid on;box on
m_contour(lg_bat,lt_bat,z,[-4000:500:0],'color',[.5 .5 .5])
m_quiver([lon_MR(2)-2 lg_all(ind)],[lat_MR(1)+0.1 lt_all(ind)],[.2 abs(nanmean(glid_pv.v_tot_cor(1:100,ind))).*glid_pv.dir_flow(1,ind)],[0 abs(nanmean(glid_pv.v_tot_cor(1:100,ind))).*glid_pv.dir_flow(2,ind)]*cosd(nanmean(lat_MR)),5,'color',[.2 .2 .2]*0,'ShowArrowHead','off')
m_scatter(lg_all(ind),lt_all(ind),5,nanmean(SA_all(1:20,ind)),'filled')
m_text(lon_MR(2)-2,lat_MR(1)+0.3,'0.2 m s^{-1}','fontsize',5)

caxis([34.95 35.35])
if l==4
    h=colorbarnew('v',0.02,1,'S_A^{0-20} (g kg^{-1})');
end
end
colormap(ax2,cmocean('haline'))

printHR('fig/MRall_map')

% section (Vg,T,S,sigma) with all data
xl = [datenum(2016,6,1) datenum(2017,8,1)];
lx = 0.8;ly = 0.21;
g_lim = [27.5 28.06];
%
ax0=axes('position',[.1 .1+3*(ly+0.01) lx ly],'xticklabel',[]);
hold on;grid on;box on
pcolor(nanmean(DAYS_all),nanmean(glid_pv.pp,2),glid_pv.v_tot_cor);shading flat
plot(tt_all,MLD_all,'-.k','markersize',3)

zdn; ylabel('Depth (m)'); xlim(xl);ylim([0 1000])
h=colorbarnew('v',0.01,1,'V_g^{c} (m s^{-1})'); colormap(ax0,cmocean('delta')); caxis([-1 1])


ax1=axes('position',[.1 .1+2*(ly+0.01) lx ly],'xticklabel',[]);
hold on;grid on;box on
pcolor(tt_all,0:1010,CT_all);shading flat
plot(tt_all,MLD_all,'-.k','markersize',3)

zdn;ylabel('Depth (m)');xlim(xl);ylim([0 1000])
h=colorbarnew('v',0.01,1,'\Theta (\circC)');colormap(ax1,cmocean('thermal'));



ax2=axes('position',[.1 .1+1*(ly+0.01) lx ly],'xticklabel',[]);
hold on;grid on;box on
pcolor(tt_all,0:1010,SA_all);shading flat
plot(tt_all,MLD_all,'-.k','markersize',3)

zdn;ylabel('Depth (m)');xlim(xl);ylim([0 1000])
h=colorbarnew('v',0.01,1,'S_A (g kg^{-1})');colormap(ax2,cmocean('haline'))


ax3=axes('position',[.1 .1 lx ly]);
hold on;grid on;box on
pcolor(tt_all,0:1010,SIG_all);shading flat
plot(tt_all,MLD_all,'-.k','markersize',3)

zdn; ylabel('Depth (m)'); xlabel('Time');xlim(xl);ylim([0 1000])
h=colorbarnew('v',0.01,1,'\sigma_0 (kg m^{-3})');colormap(ax2,cmocean('haline'))

datetick('x','keeplimits')

axes('position',[.1 .1+4*ly+0.025 lx 0.015],'xticklabel',[]);
sec_beg = [40 153 232 299 422 445 553 707 802 930 964 1010 1057 1486 1549 1616 1678 1756 1860 1963 2021 2129 2201 2303 2366];
sec_end = [130 215 277 352 441 538 613 780 872 964 1010 1045 1100 1537 1597 1662 1756 1847 1901 2021 2056 2173 2270 2366 2428];
hold on;
for l=1:length(sec_end)
    plot(nanmean(tt_all([sec_beg(l) sec_end(l)])),0*ones(2,1),'-sb','markersize',2,'markerfacecolor','b')
    text(nanmean(tt_all([sec_beg(l) sec_end(l)]))-2,1,['\bfS' num2str(l)],'fontsize',8)
end
axis off;xlim(xl)

printHR('fig/MRall_section')

sec_beg = [40 153 232 299 422 445 553 707 802 930 964 1010 1057 1486 1549 1616 1678 1756 1860 1963 2021 2129 2201 2303 2366];
sec_end = [130 215 277 352 441 538 613 780 872 964 1010 1045 1100 1537 1597 1662 1756 1847 1901 2021 2056 2173 2270 2366 2428];
rev = [0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
% Build mask from TS space of different water masses, usefull to compute Atlantic Water transport
[sASW_all wASW_all AW_all mAW_all nAW_all NSIW_all NSDW_all] = select_water_mass(SA_all,CT_all);

% get along-track distance and along-stream distance (perpendicular to DAC)
d_along = nanmean(glid_pv.d_along); dc_along = glid_pv.d_cor;
reload_io = 0; % set to 1 to force reloading optimal interpolation
zz = dir('data');
clear glid_MR TR
if ~isempty(find(strcmp({zz.name},'glid_MR.mat')))
    load('data/glid_MR.mat')
else
    reload_io = 1; % if no matfile then force to reload optimal interpolation
end

for ll= 1:length(sec_beg)
    opt_vis = 0;
    glid_MR.Lsmo = 20;glid_MR.Zsmo = 15;glid_MR.err = 0.05; % correlation length scale for OI
    dr = 2.5;dz = 5; % grid resolution for optimal ineterpolation 
    comp_glid_section_MR % do the MR front analysis for an individual section
    plot_glid_section_MR % plot some figures % skip to process faster
end
% save results to avoid recomputing everytime
save('data/glid_MR.mat','glid_MR','TR');

glid_MR

glid_MR

% synthesis transport
for transp = {'','bt','bc'} % mass heat salt
    for varr = {'','H','S'} % mass heat salt
    % different water masses
        for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
        eval([transp{1} varr{1} 'Tr' watermass{1} ' = -[TR.' transp{1} varr{1} 'Tr' watermass{1} '];'])
        end
    end
end
month_front = str2num(datestr(tt_front,'mm'));
front_width = (glid_MR.iiq2-glid_MR.iiq1)*dr;
front_width_gauss = glid_MR.fit_gauss(:,3);

iqJFMA_all = find(month_front==1 | month_front==2 | month_front==3 | month_front==4);
iqMJJA_all = find(month_front==5 | month_front==6 | month_front==7 | month_front==8);
iqSOND_all = find(month_front==9 | month_front==10 | month_front==11 | month_front==12);

load data/zz_smo.mat
lgg = linspace(lon_MR(2),lon_MR(1),500);ltt = linspace(lat_MR(1),lat_MR(2),500);
[xx yy] = earthcoord('flat',lgg,ltt,nanmean(lg_front),nanmean(lt_front));
[xx_p yy_p] = project_pt(xx,yy,-nanmean(dir_front(1,:))/nanmean(dir_front(2,:)),0);
[lggg lttt] = xy_to_lglt(xx_p,yy_p,nanmean(lg_front),nanmean(lt_front));
batfront = interp2(lg_bat,lt_bat,z_smo,lg_front,lt_front);

col_paired = paired();

[x_front y_front] = earthcoord('flat',lg_front,lt_front,nanmean(lg_front),nanmean(lt_front));
[xfront_p yfront_p] = project_pt(x_front,y_front,-nanmean(dir_front(1,:))/nanmean(dir_front(2,:)),0);
[lgfront_p ltfront_p] = xy_to_lglt(xfront_p,yfront_p,nanmean(lg_front),nanmean(lt_front));

batsmo_mean = interp2(lg_bat,lt_bat,z_smo,lggg,lttt);bat_mean = interp2(lg_bat,lt_bat,z,lggg,lttt);
[xxx yyy] = earthcoord('flat',lggg,lttt,nanmean(lgfront_p),nanmean(ltfront_p));
ddd = [0 cumsum(sqrt(diff(xxx).^2+diff(yyy).^2))];
ddd = ddd - interp1(lttt,ddd,nanmean(ltfront_p));

% compute average section
compute_avMR

% make_MRclimato
% plot mean section
xx=dd_grid;yy=pp_grid;tt=CT_mean;ss=SA_mean;vg=VG_mean;sig=SIG_mean;vger=VG_err;opt='_zmean';nnobs=nobs;plot_meanMR
for k=1:3
tt=ct_seas(:,:,k);ss=sa_seas(:,:,k);vg=vg_seas(:,:,k);sig=sig_seas(:,:,k);vger=vg_err_seas(:,:,k);opt=['_zmean_seas' num2str(k)];nnobs=nobs_seas(:,:,k);plot_meanMR
end

% figure 1 similar of Bosse and Fer, GRL 2019
plot_paper_MRmap

% figure 2 similar of Bosse and Fer, GRL 2019 and other figures (seasonal front strcuture and water mass detection)
plot_paper_MRmeansection

% display some front characteristics and build latex table 
get_PF_charac

% figure 3
plot_paper_MRtseries
