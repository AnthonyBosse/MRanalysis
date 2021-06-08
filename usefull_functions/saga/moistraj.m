function [M] = moistraj(debut,fin,suf,clip,rep1,rep2,b)
n  = length(suf);
%
for x=1:n
	load([rep1 'fpos' suf(x,:)])
	load([rep2 'fjptw' suf(x,:)])
end
for x=1:n;
	ptw = eval(['fjptw' suf(x,:)]);
	ptw(:,3) = filtfilt(b,1,ptw(:,3));
	ptw(:,2) = filtfilt(b,1,ptw(:,2));
	for i=1:size(ptw(:,3)),
		tp(i) = sw_ptmp(38.5,ptw(i,3),ptw(i,2),10);
	end;
	po = eval(['fpos' suf(x,:)]);
	caxis(clip)
	scatter(po(debut:fin,2),po(debut:fin,3),8,tp(debut:fin));
	hold on
end
%______________________________
%
%	cotes - legende
%______________________________
%
colorbar('vert')
cotmed(39.5,43.5,3,7,'k',1.5);
grid on
hold on
xlabel('longitude EST')
ylabel('latitude NORD')

