function [] = plot_m_triangle(x,y,dim,col)

tmp = get(gcf,'papersize');

xx = [0 -sind(60)*dim sind(60)*dim 0];
yy = [dim -cosd(60)*dim -cosd(60)*dim dim];

theta = 90 -  (double((x(end)-x(end-1))<0)*180 + atan( (y(end)-y(end-1))/cosd(y(end)) / (x(end)-x(end-1)) )*180/pi);

M = [cosd(theta) sind(theta) ; -sind(theta) cosd(theta)];
for l=1:length(xx)
tmp = M*[xx(l); yy(l)];
xxx(l) = tmp(1);yyy(l) = tmp(2);
end

m_patch(xxx+x(end),yyy*cosd(y(end))+y(end),col,'linestyle','none')
