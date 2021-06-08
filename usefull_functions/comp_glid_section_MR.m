% define front indices, flip left-right the reverse sections and extract relevant variables
ind_front = sec_beg(ll):sec_end(ll);

if rev(ll) == 0
    d_sec = d_along(ind_front)-d_along(ind_front(1));dd_sec = dc_along(ind_front)-dc_along(ind_front(1));
    v_tot_tmp = glid_pv.v_tot_cor(:,ind_front); % with along-stream geometric correction
    %v_tot_tmp = glid_pv.v_tot(:,ind_front);dd_sec = d_sec; % keep cross-track distance
else
    d_sec = -fliplr(d_along(ind_front)-d_along(ind_front(end)));dd_sec = -fliplr(dc_along(ind_front)-dc_along(ind_front(end)));
    ind_front = fliplr(ind_front);
    v_tot_tmp = -glid_pv.v_tot_cor(:,ind_front); % with along-stream geometric correction
    %v_tot_tmp = -glid_pv.v_tot(:,ind_front);dd_sec = d_sec; % keep cross-track distance
end
ct_tmp = CT_all(:,ind_front);pp_tmp = PP_all(:,ind_front);
sa_tmp = SA_all(:,ind_front);ssa_tmp = glid_pv.SAs(:,ind_front);
sig_tmp = SIG_all(:,ind_front);sigs_tmp = glid_pv.SIGs(:,ind_front);

% pv and insta
pvA_tmp = glid_pv.pv_A_cor(:,ind_front);pvB_tmp = glid_pv.pv_B_cor(:,ind_front);
pv_tmp = glid_pv.pv_cor(:,ind_front);
M4_tmp = (glid_pv.dBdx_cor(:,ind_front)).^2;N2_tmp = glid_pv.N2(:,ind_front);
vort_tmp = glid_pv.vort_cor(:,ind_front);ff_tmp = glid_pv.ff(:,ind_front);

theta_tmp = nanmean(glid_pv.sint(:,ind_front));
llg_glid = lg_all(ind_front); llt_glid = lt_all(ind_front); tt_glid = tt_all(ind_front);
ddir_front = glid_pv.dir_flow(:,ind_front);
for watermass = {'sASW','wASW','AW','nAW','mAW','NSIW','NSDW'}
    eval([watermass{1} '_tmp = ' watermass{1} '_all(:,ind_front);'])
    eval([watermass{1} '_tmp = ' watermass{1} '_tmp(1:5:end-6,:);'])
end

ct_tmpv = ct_tmp(1:5:end-6,:);sa_tmpv = sa_tmp(1:5:end-6,:);sig_tmpv = sig_tmp(1:5:end-6,:);


% detect the front jet in profiles of geostrophic velocities
vsurf = nanmean(v_tot_tmp(1:10,:)); % average velocity in 0-50m
if ll==22 % remove first part of data, eddy...
    vsurf(1:15) = NaN;
end
[vm iq_front] = nanmin(vsurf); % iq_front = where vsurf < 0 around peak velocity

% get front indices where v < -0.02m/s around peak velocity
indneg = find(vsurf<-0.02);
dindneg = diff(indneg);
l = find(indneg == iq_front);
front_ind = [l];
while l<length(dindneg) && dindneg(l) == 1
    front_ind = [front_ind l+1];
    l = l+1;
end
l = find(indneg == iq_front)-1;
while l>0 && dindneg(l) == 1
    front_ind = [l front_ind];
    l = l-1;
end
i_front(ll) = iq_front;
d_sec = d_sec - d_sec(iq_front);dd_sec = dd_sec - dd_sec(iq_front);

% some front characteristics (direction and position)
lg_front(ll) = llg_glid(iq_front); lt_front(ll) = llt_glid(iq_front); tt_front(ll) = tt_glid(iq_front);
front_inddir = find(abs(d_sec(indneg))<10);
dir_front(:,ll) = wmean(ddir_front(:,indneg(front_inddir))',abs([nanmean(v_tot_tmp(1:20,indneg(front_inddir)))' nanmean(v_tot_tmp(1:20,indneg(front_inddir)))']));


% Step 1 : bin CT,SA in cross-front axis (grid defined by dr and dz)
% build meshgrid
dpro_sca = repmat(dd_sec,size(ct_tmp,1),1);
dv_sca = repmat(dd_sec,size(v_tot_tmp,1),1);
bin_d = -100:dr:100; % vec bin en km
bin_p = 0:dz:1000;
[d_grid p_grid] = meshgrid(bin_d,bin_p);
dd_grid = red_xy(d_grid);pp_grid = red_xy(p_grid);
glid_MR.bin_d = dd_grid;glid_MR.bin_p = pp_grid;
% start binning and optimal interpolation
vv_bin = bin1dr(dd_sec,nanmean(v_tot_tmp),glid_MR.bin_d(1,:));
glid_MR.vv_bin(ll,:) = vv_bin.mean;
    
if reload_io
    disp('Binning TS')
    ct_bin = bin2dr(dpro_sca,pp_tmp,ct_tmp,dd_grid,pp_grid,0);
    sa_bin = bin2dr(dpro_sca,pp_tmp,sa_tmp,dd_grid,pp_grid,0);
    disp('Optimal interpolation TS')
    iqnan = find(~isnan(ct_bin.mean));
    [ct_io,ct_err]=objmap(dd_grid(iqnan),pp_grid(iqnan),ct_bin.mean(iqnan),dd_grid,pp_grid,[glid_MR.Lsmo glid_MR.Zsmo],glid_MR.err); % with default parameters and verbose = off
    ct_io(ct_err>2*glid_MR.err) = NaN;glid_MR.CT_io(:,:,ll) = ct_io;glid_MR.CT_err(:,:,ll) = ct_err;
    iqnan = find(~isnan(sa_bin.mean));
    [sa_io,sa_err]=objmap(dd_grid(iqnan),pp_grid(iqnan),sa_bin.mean(iqnan),dd_grid,pp_grid,[glid_MR.Lsmo glid_MR.Zsmo],glid_MR.err);
    sa_io(sa_err>2*glid_MR.err) = NaN;glid_MR.SA_io(:,:,ll) = sa_io;glid_MR.SA_err(:,:,ll) = sa_err;
    
    % Step 2 : reconstruct Vgeo from optimal inteprolation of CT,SA reference with DAC
    geo_strf = gsw_geo_strf_dyn_height(glid_MR.SA_io(:,:,ll),glid_MR.CT_io(:,:,ll),glid_MR.bin_p,0);
    [v_geo XXX] = gradient(geo_strf/gsw_f(lt_front(ll)),dr*1000,dz);
    glid_MR.VG_io(:,:,ll) = v_geo - repmat(nanmean(v_geo) - naninterp(glid_MR.vv_bin(ll,:)')',size(glid_MR.bin_d,1),1);
    % get deep component of the flow 
    indz = find(pp_grid(:,1)>950);
    glid_MR.VGbt_io(:,:,ll) = repmat(nanmean(glid_MR.VG_io(indz,:,ll)),size(glid_MR.bin_d,1),1); % abs vel at 1000m
    glid_MR.VGbc_io(:,:,ll) = glid_MR.VG_io(:,:,ll)-glid_MR.VGbt_io(:,:,ll);
end
glid_MR.SIG_io(:,:,ll) = gsw_sigma0(glid_MR.SA_io(:,:,ll),glid_MR.CT_io(:,:,ll));


% Step 3 : Get Front diag and transport from optimal interplation
vsurf_io = nanmean(glid_MR.VG_io(1:25,:,ll));
[vm_io iq_front_io] = nanmin(vsurf_io);
indneg_io = find(vsurf_io<-.02);
dindneg_io = diff(indneg_io);
l = find(indneg_io == iq_front_io);
front_ind_io = [l];
while l<length(dindneg_io) && dindneg_io(l) == 1
    front_ind_io = [front_ind_io l+1];
    l = l+1;
end
l = find(indneg_io == iq_front_io)-1;
while l>0 && dindneg_io(l) == 1
    front_ind_io = [l front_ind_io];
    l = l-1;
end

% Step 4 : Compute Volume/heat/salt transport for total/geostrophic shear/deep component and for different water masses   
[mask_sASW mask_wASW mask_AW mask_mAW mask_nAW mask_NSIW mask_NSDW] = select_water_mass(glid_MR.SA_io(:,:,ll),glid_MR.CT_io(:,:,ll));
Tref=-.1;Sref=35;tmpvar=ones(size(glid_MR.CT_io(:,:,ll)))/1e6;tmpvarH=(glid_MR.CT_io(:,:,ll)-Tref)*1028*gsw_cp0;tmpvarS=(Sref-glid_MR.SA_io(:,:,ll))*1028/1000/Sref;mask_ = ones(size(mask_AW));

for transp = {'','bc','bt'} % tot baroclinic bartropic with limits set frontal detection (V<-0.02m/s) (TR.*_ideal)
    for varr = {'','H','S'} % mass heat salt
        for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}%%% different water masses

        eval(['TR.H_' watermass{1} '(ll,:) = sum(mask_' watermass{1} ')*dz;'])
        eval(['tmp=trapz(glid_MR.bin_p(:,1),naninterp2_extrap(glid_MR.VG' transp{1} '_io(:,:,ll).*mask_' watermass{1} '.*tmpvar' varr{1} '));tmp(isnan(tmp)) = 0;' transp{1} varr{1} 'Tr' watermass{1} '_int =cumtrapz(glid_MR.bin_d(1,:)*1000,tmp);'])
        eval(['TR.' transp{1} varr{1} 'Tr' watermass{1} '_int(ll,:) = ' transp{1} varr{1} 'Tr' watermass{1}  '_int - ' transp{1} varr{1} 'Tr' watermass{1}  '_int(indneg_io(front_ind_io(1))); TR.' transp{1} varr{1} 'Tr' watermass{1}  '(ll) = TR.' transp{1} varr{1} 'Tr' watermass{1}  '_int(ll,indneg_io(front_ind_io(end)));'])

        end
    end
end

% some front characteristics (width, intensity, ...) and store into a strcuture glid_MR
glid_MR.front_width(ll) = glid_MR.bin_d(1,indneg_io(front_ind_io(end))) - glid_MR.bin_d(1,indneg_io(front_ind_io(1))); 
glid_MR.front_intens(ll) = abs(vm_io);
glid_MR.vsurf{ll} = vsurf_io;
glid_MR.indneg{ll} = indneg_io;
glid_MR.front_ind{ll} = front_ind_io;
glid_MR.iq_front{ll} = iq_front_io;
glid_MR.lg_front(ll) = llg_glid(iq_front);
glid_MR.lt_front(ll) = llt_glid(iq_front);
glid_MR.tt_front(ll) = tt_glid(iq_front);
for watermass = {'AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
    eval(['glid_MR.mask_' watermass{1} '{ll} = mask_' watermass{1} ';'])
end
glid_MR.dir_front(:,ll) = wmean(ddir_front(:,indneg(front_inddir))',abs([nanmean(v_tot_tmp(1:20,indneg(front_inddir)))' nanmean(v_tot_tmp(1:20,indneg(front_inddir)))']));	


%%%% Ste 5: Gaussian fitting
%%%% fit a gaussian jet, deep bt
Vjet = nanmean(glid_MR.VG_io(180:end,:,ll)); Djet = red_x(bin_d);
iqtest = abs(Vjet./nanmin(Vjet(abs(Djet)<20)) )>0.5 & abs(Djet)<20;
global Djet Vjet iqtest
x_obs = fminsearch('cost_gaussian_MRjet',[0 vm_io 20]);
x_obs(2) = nanmin(Vjet(iqtest));
Vfitgbt = x_obs(2)*gauss_distribution(Djet, x_obs(1), x_obs(3));

%%%% Fit a gaussian jet on average 0-1000m velocity
Vjet = nanmean(glid_MR.VG_io(:,:,ll)); Djet = red_x(bin_d);
iqtest = abs(Vjet./nanmin(Vjet) )>0.5 & abs(Djet)<20;
global Djet Vjet iqtest
x_obs = fminsearch('cost_gaussian_MRjet',[0 vm_io 20]);
x_obs(2) = nanmin(Vjet(iqtest));
Vfitg = x_obs(2)*gauss_distribution(Djet, x_obs(1), x_obs(3));


% deal with particular cases
if ll==20
    [v iq1]=nanmin(abs(Djet-(x_obs(1)-x_obs(3))));[v iq2]=nanmin(abs(Djet-(x_obs(1)+0.5*x_obs(3))));
elseif ll==23
    [v iq1]=nanmin(abs(Djet-(x_obs(1)-1*x_obs(3))));[v iq2]=nanmin(abs(Djet-(x_obs(1)+2*x_obs(3))));
elseif ll==9
    [v iq1]=nanmin(abs(Djet-(x_obs(1)-2*x_obs(3))));[v iq2]=nanmin(abs(Djet-(x_obs(1)+x_obs(3))));
else
    [v iq1]=nanmin(abs(Djet-(x_obs(1)-2*x_obs(3))));[v iq2]=nanmin(abs(Djet-(x_obs(1)+2*x_obs(3))));
end
% iq1/2 -> width given by guassian fit with +/- 1 to 2 sigma depending on sections
% then take the min of theoretical width and width of front (<0 vel) deduced from optimal interppolation
iiq1 = max(iq1,indneg_io(front_ind_io(1)));iiq2 = min(iq2,indneg_io(front_ind_io(end)));

glid_MR.iiq1(ll) = iiq1;glid_MR.iiq2(ll) = iiq2;
glid_MR.fit_gauss(ll,:) = x_obs;

disp(['## Front ' num2str(ll) '/' num2str(length(sec_beg)) ' : optimal interpolation Vmax = ' num2str(vm_io,3) ' m/s, width ' num2str(glid_MR.front_width(ll),2) ' km , DAC Gaussian fit (' num2str(x_obs(2),2) ' m/s, ' num2str(x_obs(3),2) ' km)'])

% Volume transport with limits set by gaussian fit (TR.*_int_gauss)
for transp = {'','bc','bt'} % tot baroclinic bartropic
    for varr = {'','H','S'} % mass heat salt
        %%% different water masses
        for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
        eval(['TR.' transp{1} varr{1} 'Tr' watermass{1} '_int_gauss(ll,:) = ' transp{1} varr{1} 'Tr' watermass{1}  '_int - ' transp{1} varr{1} 'Tr' watermass{1}  '_int(iiq1); TR.' transp{1} varr{1} 'Tr' watermass{1}  '(ll) = TR.' transp{1} varr{1} 'Tr' watermass{1}  '_int_gauss(ll,iiq2);'])
        end
    end
end

%%% transport with limits set by gaussian fit and idealized velocities (TR.*_ideal)
for transp = {'','bt'}
    for varr = {'','H','S'} % mass heat salt
        for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
        eval(['TR.' transp{1} varr{1} 'Tr' watermass{1} '_ideal(ll) = trapz(glid_MR.bin_d(1,iiq1:iiq2)*1000,TR.H_' watermass{1} '(ll,iiq1:iiq2).*Vfitg' transp{1} '(iiq1:iiq2).*nanmean(tmpvar' varr{1} '(:,iiq1:iiq2)));'])
        end
    end
end
for varr = {'','H','S'} % mass heat salt
    for watermass = {'','AW','mAW','nAW','NSDW','NSIW','sASW','wASW'}
    eval(['TR.bc' varr{1} 'Tr' watermass{1} '_ideal(ll) = TR.Tr' watermass{1} '_ideal(ll) - TR.btTr' watermass{1} '_ideal(ll);'])
    end
end

Vmax_gauss(ll) = x_obs(2);
Vsurf_gauss(ll) = nanmin(nanmin(glid_MR.VG_io(1:5,iqtest,ll)));
Vbot_gauss(ll) = nanmin(nanmin(glid_MR.VG_io(end-5:end,iqtest,ll)));
front_width_gauss(ll) = x_obs(3);


dt_sec(ll) = diff(mima(tt_glid));
dt_front(ll) = diff(mima(tt_glid(front_ind)));

disp(['### mean transport from optimal interpolation : mean ' num2str(TR.Tr(ll),3) ' Sv -- mean barotropic ' num2str(TR.btTr(ll)) ' Sv -- mean baroclinic ' num2str(TR.bcTr(ll),3) ' -- Tr AW mean ' num2str(TR.TrAW(ll),3)])
if 0
%%%% error estimates
for krand = 1:100
[mask_sASW mask_wASW mask_AW mask_mAW mask_nAW mask_NSIW mask_NSDW] = select_water_mass(glid_MR.SA_io(:,:,ll),glid_MR.CT_io(:,:,ll));mask_ASW = (mask_wASW+mask_sASW)/2;
Tref=-.1;Sref=35;tmpvar=ones(size(glid_MR.CT_io(:,:,ll)))/1e6;tmpvarH=(glid_MR.CT_io(:,:,ll)-Tref)*1028*gsw_cp0;tmpvarS=(Sref-glid_MR.SA_io(:,:,ll))*1028/1000/Sref;mask_ = ones(size(mask_AW));
U_err = 0.02*repmat(rand(1,size(mask_AW,2))*2-1,size(mask_AW,1),1)+(rand(size(mask_AW))*2-1).*glid_MR.VG_err(:,:,ll);U_errbc = 0.02*repmat(rand(1,size(mask_AW,2))*2-1,size(mask_AW,1),1);
for transp = {'','bc','bt'} % tot baroclinic bartropic
for varr = {'','H'} % mass heat salt
for watermass = {'','AW','mAW','nAW','ASW'}%%% different water masses
if strcmp(transp(1),'bc')
U_err = U_err + U_errbc;
elseif strcmp(transp(1),'bt')
U_err = U_err - U_errbc;
end
eval(['tmp=trapz(glid_MR.bin_p(:,1),naninterp2_extrap((glid_MR.VG' transp{1} '_io(:,:,ll)+U_err).*mask_' watermass{1} '.*tmpvar' varr{1} '));tmp(isnan(tmp)) = 0;' transp{1} varr{1} 'Tr' watermass{1} '_int =cumtrapz(glid_MR.bin_d(1,:)*1000,tmp);'])

eval(['ERR.' transp{1} varr{1} 'Tr' watermass{1} '_int(ll,:) = ' transp{1} varr{1} 'Tr' watermass{1}  '_int - ' transp{1} varr{1} 'Tr' watermass{1}  '_int(iiq1); ERR.' transp{1} varr{1} 'Tr' watermass{1}  '(ll,krand) = ERR.' transp{1} varr{1} 'Tr' watermass{1}  '_int(ll,iiq2);'])
end
end
end
end
disp(['### stat error (n=' num2str(krand) '), error 2 cm/s on DAC : mean ' num2str(nanmean(ERR.Tr(ll,:)),3) ' std ' num2str(nanstd(ERR.Tr(ll,:)),3) ' -- mean barotropic ' num2str(nanmean(ERR.btTr(ll,:)),3) ' std ' num2str(nanstd(ERR.btTr(ll,:)),3) ' -- mean baroclinic ' num2str(nanmean(ERR.bcTr(ll,:)),3) ' std ' num2str(nanstd(ERR.bcTr(ll,:)),3) ' -- Tr AW mean ' num2str(nanmean(ERR.TrAW(ll,:)),3) ' std ' num2str(nanstd(ERR.TrAW(ll,:)),3)])
stat_Tr(ll) = nanmean(ERR.Tr(ll,:));stat_Tr_er(ll) = nanstd(ERR.Tr(ll,:));
stat_TrAW(ll) = nanmean(ERR.TrAW(ll,:));stat_TrAW_er(ll) = nanstd(ERR.TrAW(ll,:));
end
