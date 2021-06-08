function [x3,y3]=earthcoord(proj,lg,lt,lgc,ltc)
%[x3,y3]=earthcoord(proj,lg,lt,lgc,ltc)
% proj = 'ster'
%     function [x3,y3]=coord_ster(lg,lt,lgc,ltc,vect)
%     x3 et y3 : matrice (grille) des coordonnees
%                stereographiques polaires
%                taille [length(x1),length(y1)]
%     vect=parametre optionnel : 1=> x3 et y3=vecteurs (x1,y1 de meme taille)
%     Origine en (lgcE,90N)
%     x1=longitude, y1=latitude (vecteurs ou matrices)
%     y1 : pour l'hemisphere Sud, utiliser coordster(lg,-lt,lgc).
%
% proj = 'merc'
%     function [x3,y3]=coordmerc(lg,lt,lgc,ltc,vect)
%     x3 et y3 : matrice (grille) des coordonnees 
%           en projection de Mercator (aspect direct)
%           taille [length(x1),length(y1)]
%     vect : parametre optionnel. 1 => x3 et y3 sont des vecteurs
%     x1=longitude, y1=latitude (vecteurs ou matrices)
%     y1 : pour l'hemisphere Sud, utiliser coordmerc(lg,-lt).
% proj = 'flat'
%     projection plate en km 
%
%
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
  % todo : scale in km for mercator / (lgc,ltc)
end


