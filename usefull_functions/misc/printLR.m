function [] = printLR(filename)

v = version('-release');
if str2num(v(1:4))<2014
eval(['print -depsc2 ' filename '.eps'])
eval(['!epstopdf ' filename '.eps'])
eval(['!rm ' filename '.eps'])
eval(['!convert -density 300 ' filename '.pdf ' filename '.png'])
else
set(gcf,'color','w');
def_fontname
eval(['export_fig -r200 ' filename '.png'])
%eval(['export_fig ' filename '.pdf'])
%print(['print -dpng -r200 ' filename '.png'])
end
