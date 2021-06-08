function dmy=julien(jd);

%
%	dmy = julien (jd)
%
%	convertir des jours juliens
%	dans le calendrier gregorien
%
%	retourne matrice de trois colonnes (dmy)
%	contenant le jours, le mois et l'annee
%	correspondant au jours julien (jd)
%


%
%	permet de convertir les jd utilises pour
%	sofargos (9xxx) et pour esop (xxx)
%

jd=jd(:);
nj=find(jd<4000);	% cad jj sup a 2450000 3 digit
jd(nj)=jd(nj)+2450000.5;
nj=find(jd<10000);	% cad jj sup a 2440000 4 digit
jd(nj)=jd(nj)+2440000.5;

z=fix(jd+0.5);
f=jd+.5-z;
al=fix((z-1867216.25)/36524.25);
a=z+1+al-fix(al/4);
nz=find(z<2299161);
a(nz)=z(nz);
b=a+1524;
c=fix((b-122.1)/365.25);
d=fix(365.25*c);
e=fix((b-d)/30.6001);
d=b-d-fix(30.6001*e)+f;
m=e-1;
ne=find(e>13);
m(ne)=e(ne)-13;
y=c-4715;
nm=find(m>2);
y(nm)=c(nm)-4716;
dmy=[d m y];
