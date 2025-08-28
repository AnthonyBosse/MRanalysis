function [glid_pv] = comp_runningPV(glid_stru,Lsmo,dx)
% input : glid_stru = glider structure with following fields
%         Lsmo = filtering length scale Lsmo (gaussian variance)
%         dx = 10; % if dx>0, grid the data on a regular grid (bin average), otherwise keep original profiles
%
% Example of required data structure:
%        glider_name: 'sg560'
%         deployment: 'MR1'
% Vertical profiles : pressure PP, time DAYS, position lon-lat ZLG-ZLT, cons. temperature CT, abs. salinity SA, pot. density SIG (from Gibbs Seawater toolbox)
%                 PP: [1011x1340 double]
%               DAYS: [1011x1340 double]
%                ZLG: [1011x1340 double]
%                ZLT: [1011x1340 double]
%                 CT: [1011x1340 double]
%                 SA: [1011x1340 double]
%                SIG: [1011x1340 double]
% Dive-average currents : position lon-lat LGV-LTV, zonal and meridionnal component VX-VY, time TVEL
%                LGV: [1x673 double]
%                LTV: [1x673 double]
%                 VX: [1x673 double]
%                 VY: [1x673 double]
%               TVEL: [1x673 double]
%
% Anthony Bosse, December 2019 (anthony.bosse@mio.osupytheas.fr)
% last update, November 2024 (allows output on regular grid setting dx param >0)

warning off
disp(['### compute PV for ' glid_stru.glider_name '(' glid_stru.deployment ')'])

%%% get along-track distance
lg_ref = nanmean(glid_stru.ZLG(:));lt_ref = nanmean(glid_stru.ZLT(:));
[xxkm yykm] = earthcoord('flat',nanmean(glid_stru.ZLG),nanmean(glid_stru.ZLT),lg_ref,lt_ref);
d_along = NaN_interp(interp1(nanmean(glid_stru.DAYS),[0 nancumsum(sqrt(diff(xxkm).^2 + diff(yykm).^2))],glid_stru.DAYS(:)),1);
xxkm = NaN_interp(interp1(nanmean(glid_stru.DAYS),xxkm,glid_stru.DAYS(:)),1);yykm = NaN_interp(interp1(nanmean(glid_stru.DAYS),yykm,glid_stru.DAYS(:)),1);

xx_km = reshape(xxkm,size(glid_stru.ZLG));
yy_km = reshape(yykm,size(glid_stru.ZLG));
d_along = reshape(d_along,size(glid_stru.ZLG));
d_along=d_along-d_along(1); % reset first profile to 0 km

%%% interp dive-average currents (DAC) on profiles time
dir_x = gradient(nanmean(xx_km),nanmean(d_along));
dir_y = gradient(nanmean(yy_km),nanmean(d_along));
dir_norm = sqrt(dir_x.^2+dir_y.^2);
dir_x = dir_x./dir_norm;dir_y = dir_y./dir_norm; % normalized along-track vectors

[xxkmvel yykmvel] = earthcoord('flat',glid_stru.LTV,glid_stru.LGV,lg_ref,lt_ref);

vx_pro = interp1d(glid_stru.TVEL,glid_stru.VX,nanmean(glid_stru.DAYS));vy_pro = interp1d(glid_stru.TVEL,glid_stru.VY,nanmean(glid_stru.DAYS));
vxn_pro=vx_pro./sqrt(vx_pro.^2+vy_pro.^2);vyn_pro=vy_pro./sqrt(vx_pro.^2+vy_pro.^2);

%%% angle between glider track and DAC
theta = real(acosd(dot([dir_x ; dir_y],[vxn_pro ; vyn_pro])));
theta_smo = mean_av_gauss(theta,nanmean(d_along),Lsmo);
theta(theta==0)=NaN;%theta = naninterp(theta);
sintheta = sind(theta_smo);

% QC on angle
indgood = find(abs(sintheta)>0.2);%indgood = 1:length(vx_pro); % optionnal, to exclude when angle is too large, then interp over excluded values
vx_pro = interp1d(nanmean(glid_stru.DAYS(:,indgood)),vx_pro(indgood),nanmean(glid_stru.DAYS));
vy_pro = interp1d(nanmean(glid_stru.DAYS(:,indgood)),vy_pro(indgood),nanmean(glid_stru.DAYS));
sintheta = interp1(indgood,sintheta(indgood),1:length(sintheta));

%%%% Gridding and smoothing of CT,SA sections before computing gradients
dz = 5; % interp on vertical every dz meters before computing diags
if dx==0
[bin_d bin_p] = meshgrid(nanmean(d_along),(0:dz:1000)');
SAb = griddata(d_along(:,indgood),repmat((1:size(glid_stru.PP,1))',1,length(indgood)),glid_stru.SA(:,indgood),bin_d,bin_p);
CTb = griddata(d_along(:,indgood),repmat((1:size(glid_stru.PP,1))',1,length(indgood)),glid_stru.CT(:,indgood),bin_d,bin_p);
SIGb = gsw_sigma0(SAb,CTb);
else % requires bin2d package for bin-averaging
[bin_d bin_p] = meshgrid(0:dx:nanmax(d_along(:)),(0:dz:1000)');
CTbin = bin2d(d_along(:,indgood),repmat((1:size(glid_stru.PP,1))',1,length(indgood)),glid_stru.CT(:,indgood),bin_d,bin_p);CTb = CTbin.mean;
SAbin = bin2d(d_along(:,indgood),repmat((1:size(glid_stru.PP,1))',1,length(indgood)),glid_stru.SA(:,indgood),bin_d,bin_p);SAb = SAbin.mean;
SIGb = gsw_sigma0(SAb,CTb);
end
disp(['smoothing hor/vert : gaussian ' num2str(Lsmo) 'km / 5m'])
SIGs = mean_av_gauss2(SIGb,nanmean(bin_d),Lsmo);
BBsmo = -9.81/1028*(1000+SIGs);iqnan = isnan(BBsmo); % hor smooth
Bzsmo = naninterp2(mean_av_gauss2(BBsmo',bin_p(:,1)',5)');Bzsmo(iqnan) = NaN; % vert smooth
[dBdx N2] = gradient(Bzsmo,nanmean(bin_d)*1000,-nanmean(bin_p,2));

ff = gsw_f(repmat(nanmean(glid_stru.ZLT(:)),size(bin_p,1),1));
N = real(sqrt(N2));
N(size(N,1)-round(20/dz):end,:) = NaN;N(1:round(20/dz),:) = NaN;

%%%% Compute geostrophic velocities from dynamic height (from hor smoothed of T and S)
disp(['computing geostr from dyn height: smoothing gaussian ' num2str(Lsmo) 'km'])
SAs = mean_av_gauss2(SAb,nanmean(bin_d),Lsmo);CTs = mean_av_gauss2(CTb,nanmean(bin_d),Lsmo);
geo_strf = gsw_geo_strf_dyn_height(SAs,CTs,bin_p,0);
[v_geo XXX] = gradient(geo_strf./ff,nanmean(bin_d)*1000,-nanmean(bin_p,2));

%%%%% get ref velocity from DAC
v_ref = dot([-dir_y ; dir_x],[vx_pro ; vy_pro]);
v_ref = mean_av_gauss(v_ref,nanmean(d_along),Lsmo);
if dx>0 % if dx>0, interpolate DAC on a regular grid
v_ref = interp1d(nanmean(d_along),v_ref,nanmean(bin_d));
end
v_tot = v_geo + repmat(v_ref-nanmean(v_geo),size(bin_p,1),1); % total geostrophic vel, orthogonal to DAC

%%%% rel vorticity and PV
[vort XXX] = gradient(v_tot,nanmean(bin_d)*1000,-nanmean(bin_p,2)); % old style vort2 = diff(DD_vg.*vcyclo,1,2)/1000./red_x(DD_vg)./nanmean(diff(bin_d));
if all(isnan(nanmean(vort))); vort = zeros(size(vort)); end % if no DAC
[XXX dVdz] = gradient(v_tot,nanmean(bin_d)*1000,-nanmean(bin_p,2));
pv_A = (ff + vort).*N2/9.81; pv_B = - dBdx.*dVdz/9.81;
if all(isnan(nanmean(pv_B))); pv_B = -dBdx.^2./ff/9.81; end % if no DAC
pv = pv_A + pv_B ;

%%%% apply geometrical correction to get along-stream (along-DAC) vel and PV
if dx>0 % if dx>0, interpolate on regular grid
sintheta = interp1d(nanmean(d_along),sintheta,nanmean(bin_d));
sintheta(1) = 1; % avoid NaN in first element due to interp
end
sint = repmat(sintheta,size(bin_p,1),1);
pv_A_cor = (ff + vort./sint./sint).*N2/9.81;pv_B_cor = - dBdx.*dVdz/9.81./sint./sint;
pv_cor = pv_A_cor + pv_B_cor;
%
v_ref_cor = sign(dot([-dir_y ; dir_x],[vx_pro ; vy_pro])).*sqrt(vx_pro.^2 + vy_pro.^2);
v_ref_cor = mean_av_gauss(v_ref_cor,nanmean(d_along),Lsmo);
if dx>0 % if dx>0, interpolate DAC on a regular grid
v_ref_cor = interp1d(nanmean(d_along),v_ref_cor,nanmean(bin_d));
end
v_tot_cor = v_geo./sint + repmat(v_ref_cor-nanmean(v_geo./sint),size(bin_p,1),1);
vort_cor = vort./sint./sint;

% store results in structure
glid_pv.glider_name = glid_stru.glider_name;
glid_pv.deployment = glid_stru.deployment;
glid_pv.dir_flow = [vxn_pro ; vyn_pro]; % direction of DAC (normalized)
dac_sign = sign(dot([-dir_y ; dir_x],[vxn_pro ; vyn_pro]));
glid_pv.dir_DACflow = [vxn_pro.*dac_sign  ; vyn_pro.*dac_sign]; % direction of cross-stream flow (normalized)
glid_pv.dir_crosstrack = [-dir_y ; dir_x]; % direction of cross-track trajectory (normalized)
glid_pv.pp = bin_p(:,1);
glid_pv.d_along = bin_d(1,:); % distance along-track
glid_pv.d_along_cor = nancumsum([0 diff(glid_pv.d_along).*red_x(sintheta)]);
glid_pv.CTs = CTs;glid_pv.SAs = SAs;glid_pv.SIGs = SIGs; % smoothed T,S,rho used in thermal wind
glid_pv.sint = sint; % geometrical correction factor (angle between DAC and track)
glid_pv.v_tot = v_tot; glid_pv.v_tot_cor = v_tot./sint; % across-track and along-DAC geotrophic vel
glid_pv.ff = ff; glid_pv.N2 = N2; % coriolis and buoyancy freq
glid_pv.dBdx = dBdx; glid_pv.dBdx_cor = dBdx./sint; % hor buoyancy gradient without and with corr
glid_pv.vort = vort; glid_pv.vort_cor = vort./sint./sint; % rel vorticity without and with corr
glid_pv.pv_A = pv_A;glid_pv.pv_B = pv_B; glid_pv.pv = pv; % PV without geometrical correction
glid_pv.pv_A_cor = pv_A_cor;glid_pv.pv_B_cor = pv_B_cor; glid_pv.pv_cor = pv_cor; % PV with correction











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% some usefull functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% deal with NaN while interpolating
function [r] = interp1d(x,y,xi)

[xu ind tmp] = unique(x);
iq = find(~isnan(xu) & ~isnan(y(ind)));
if ~isempty(iq)
r = interp1(xu(iq),y(ind(iq)),xi);
else
r = nan(size(xi));
end

%%%%% reduce by one the dimension of 2D array by averaging
function B = red_x(A)
B = (A(:,1:end-1)+A(:,2:end))/2;

%%%%% convert lg,lt into km
function [x3,y3]=earthcoord(proj,lg,lt,lgc,ltc)
   r0 = 6371; % earth radius
   l = pi/180;
if strcmp(proj,'ster'),
  if isempty(lg), x3 = []; y3 = [];
  else
    if min(size(lg)) == 1,
      lg = lg(:);
      lt = lt(:);
    end;
    signy = sign(mean(lt(~isnan(lt))));
    lt = abs(lt);   
    lg = lg-lgc;    
    r  = 2*r0*tan(pi/4-lt'*l/2);
    y2 = 2*r0*tan(pi/4-ltc*l/2);    
    x3 = sin(lg*l).*r';
    y3 = signy*(-cos(lg*l).*r'+y2);
  end;
   
elseif  strcmp(proj,'flat')

  x3 = pi*(lg - lgc)/180*r0*cos(l*ltc);
  y3 = pi*(lt - ltc)/180*r0;

  
elseif  strcmp(proj,'merc')
  
  if isempty(lg), x3 = []; y3 = [];
  else
    if min(size(lg)) == 1,
      lg = lg(:);
      lt = lt(:);
    end;
    lg = lg-lgc;
    x3  = lg*l;
    y3  = log(tan(pi/4+lt*l/2));   
  end;
end

%%%%% Gaussian function
function f = gauss_distribution(x, mu, s)
p1 = -.5 * ((x - mu)/s) .^ 2;
p2 = 1; %(s * sqrt(2*pi));
f = exp(p1) ./ p2;

%%%%% weighted average
function y = wmean(x,w,dim)

if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);

%%%%% Gaussian moving average
function [mt] = mean_av_gauss(t,d_mean,km)

for i = 1:length(t)
	if isnan(t(i))
		mt(i) = NaN;
	else
	
	clear weight_tab
	iqnan = find(abs(d_mean-d_mean(i))<5*km & ~isnan(t));
	weight_tab = gauss_distribution(d_mean(iqnan),d_mean(i),km);

        mt(i) = wmean(t(iqnan),weight_tab);
	end
end

function [mt] = mean_av_gauss2(t,d_mean,km)
% t = tableau à smoother
% d_mean position des points
% km largeur de la gaussienne (le sigma)

for i = 1:length(t(:,1))
    mt(i,:) = mean_av_gauss(t(i,:),d_mean,km);
end


function TTT = naninterp2(TT)
for k=1:length(TT(1,:))
    if length(find(~isnan(TT(:,k))))>2
           TTT(:,k) = naninterp(TT(:,k));
        %TTT(:,k) = NaN_interp(TT(:,k),1); % Ilker, used mine, to avoid linear extrap-- just extend
    else
        TTT(:,k) = TT(:,k);
    end
end

function xint=NaN_interp(x,extrapflag)
% xint=NaN_interp(x,extrapflag)
% Find the NaN values in x (not matrix), and replaces
% with the interpolated values output is xint, the same size as x.
%
% extrapflag==1 forces nearest-neighbour extrapolation for start and end
% pieces of record...
if nargin<2
    extrapflag = 0 ;
end

if all(isnan(x))
    xint=x;
   
else
   
    flag=0;
    [m,n]=size(x);
    if n>1; x=x'; flag=1; end
    x1=denan(x); y1=find(~isnan(x));
    y2=[1:length(x)]';
   
    if length(x1)<2
        xint=NaN.*x;
    else
        xint=interp1(y1,x1,y2);
        if extrapflag | strcmp(extrapflag,'extrap')
            id = isnan(xint) ;
            xint=interp1(y2(~id),xint(~id),y2,'nearest','extrap');
        end
    end
    %make sure x and xint are of same dimension
    if flag; xint=xint';end
end
