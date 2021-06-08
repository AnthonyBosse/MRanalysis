function [q] = labfig(lab,coeffx,coeffy,rev,axx)


if nargin<5
axx = gca;
end

xl = get(axx,'xlim');
yl = get(axx,'ylim');
if rev ~= 1
q=text(xl(1)-diff(xl)/10*coeffx,yl(2)+diff(yl)/10*coeffy,lab,'fontsize',9,'fontname','arial');
else
q=text(xl(1)-diff(xl)/10*coeffx,yl(1)-diff(yl)/10*coeffy,lab,'fontsize',9,'fontname','arial');
end
set(q,'parent',axx);
