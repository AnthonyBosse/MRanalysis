lx = .3; ly = .4;

close all

ax1=axes('position',[.1 .1+1*(ly+0.05) lx ly]);
hold on;grid on;box on
pcolor(xx,yy,tt);shading interp
[c,h]=contour(xx,yy,sig,[27.3:.1:27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
%try indbad = find(nnobs<5); plot(xx(indbad),yy(indbad),'.k'); end
zdn;xlim([-120 80]);
caxis([-0.5 8])
h=colorbarnew('v',0.02,1,'\Theta (\circC)');
xlabel('Cross-front distance (km)')



ax2=axes('position',[.1+1*(lx+0.15) .1+1*(ly+0.05) lx ly]);
hold on;grid on;box on
pcolor(xx,yy,ss);shading interp
[c,h]=contour(xx,yy,sig,[27.3:.1:27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
%try indbad = find(nnobs<5); plot(xx(indbad),yy(indbad),'.k'); end
zdn;xlim([-120 80]);
caxis([35.05 35.35])
h=colorbarnew('v',0.02,1,'S_A (g kg^{-1})');

ax3=axes('position',[.1+1*(lx+0.15) .1+0*(ly+0.05) lx ly]);
hold on;grid on;box on
pcolor(xx,yy,vg);shading interp
[c,h]=contour(xx,yy,sig,[27.3:.1:27.6 27.75:0.05:28.1],'k');
clabel(c,h,'fontsize',5)
%try indbad = find(vger>abs(vg)); plot(xx(indbad),yy(indbad),'.k'); end
%try indbad = find(nnobs<5); plot(xx(indbad),yy(indbad),'.k'); end
zdn;xlim([-120 80]);
caxis([-.4 .4])
h=colorbarnew('v',0.02,1,'V_g (m s^{-1})');
xlabel('Cross-front distance (km)')



colormap(ax1,cmocean('thermal'));colormap(ax2,cmocean('haline'));colormap(ax3,cmocean('delta'));

printHR(['fig/MRmeanSection' opt])