clear;
% repertoires des donnees_________________________
%
rep1 = ['/usr2/com/testor/vsofar/pos/'];
rep2 = ['/usr2/com/testor/vsofar/jptw/'];
rep3 = ['/usr2/com/testor/vsofar/hydro/'];
%_________________________________________________
%
%
%_groupes de trajectoires______________%
%_toutes trajectoires____suffixes______%
suf = ['00';'03';'04';'05';'07';'08';'09';'10';'11';'13';'14';'15';'88';'89';'90';'92';'93';'96';'97';'99'];
%_trajectoires entre 1050 et 1450 m____%
sufb = ['00';'04';'05';'88';'93'];
%_trajectoires entre 300 et 650 m______%
sufh = ['03';'07';'08';'09';'10';'11';'13';'14';'15';'89';'90';'92';'96';'97';'99'];
%_trajectoires entre 300 et 450 m______%
sufh1 = ['03';'07';'08';'10';'13';'89';'92';'96';'99'];
%_trajectoires entre 500 et 650 m______%
sufh2 = ['09';'11';'14';'15';'90';'97'];
%______________________________________%
%
n  = length(suf);
nb = length(sufb);
nh = length(sufh);
nh1 = length(sufh1);
nh2 = length(sufh2);
%
for x=1:n
	load([rep1 'fpos' suf(x,:)])
	load([rep2 'fjptw' suf(x,:)])
end
	load([rep3 'posCTD0195'])
	load([rep1 'fzposv'])
%__________________________________________________
%
% on filtre les trajectoires (passe-bas)
%__________________________________________________
%
b=fir1(64,.01);
%
 pos00 = filtfilt(b,1,fpos00);
 pos03 = filtfilt(b,1,fpos03);
 pos04 = filtfilt(b,1,fpos04); 
 pos05 = filtfilt(b,1,fpos05);
 pos07 = filtfilt(b,1,fpos07);
 pos08 = filtfilt(b,1,fpos08);
 pos09 = filtfilt(b,1,fpos09);
 pos10 = filtfilt(b,1,fpos10);
 pos11 = filtfilt(b,1,fpos11);
 pos13 = filtfilt(b,1,fpos13);
 pos14 = filtfilt(b,1,fpos14);
 pos15 = filtfilt(b,1,fpos15);
 pos88 = filtfilt(b,1,fpos88);
 pos89 = filtfilt(b,1,fpos89);
 pos90 = filtfilt(b,1,fpos90);
 pos92 = filtfilt(b,1,fpos92);
 pos93 = filtfilt(b,1,fpos93);
 pos96 = filtfilt(b,1,fpos96);
 pos97 = filtfilt(b,1,fpos97);
 pos99 = filtfilt(b,1,fpos99);
%

%
