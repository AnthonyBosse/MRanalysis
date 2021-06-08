function [phiB phiC iq_stable iq_IISI iq_SI iq_GISI iq_GI] = get_insta_from_PV(M4,ff,N2,vort)
% get the instability from PV expression
% PV = ff*(ff + vort)*N2 - M4 where M4 = (dBdx^2+dBdy^2) in geostreophic approx, M4 = f * (dBdxdVdz - dBdydUdz)

phiB = atan2d( - M4 , ff.^2.*N2 ); phiC = atan2d( - ff - vort , ff );

iq_stable = phiB>phiC;
iq_IISI = (vort<0 & phiB<phiC & phiB>-45);
iq_SI = ((vort>0 & phiB<phiC) | (vort<0 & phiB<-45)) & phiB>-90;
iq_GISI = phiB<-90 & phiB>-135;
iq_GI = phiB<-135;
