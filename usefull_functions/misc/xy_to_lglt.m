function [lgg ltt] = xy_to_lglt(xxx,yyy,lgref,ltref)


ltt = yyy/(6400*pi/180)+ltref;
lgg = xxx/(6400*pi*cos(ltref*pi/180)/180)+lgref;
