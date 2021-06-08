function [cbar] = colorbarnew(orientation,alpha,beta,tit,axx)

if nargin<5
axx = gca;
end

if strcmp(orientation,'h')
% setup the colorbar and manually set the tick labels
col = 'south';
cbar=colorbar(axx,col);
cpos=get(cbar,'Position');
cpos(3)=cpos(3)*beta; % Halve the thickness
cpos(4)=max(cpos(4)/4,0.0075); % Halve the thickness
cpos(2)=cpos(2) - alpha; % Move it down
set(cbar,'xaxislocation','top')
%cpos(1) = cpos(1)+0.05;
%cpos(3) = .2; % elongate it
set(cbar,'Position',cpos);
elseif strcmp(orientation,'v')
% setup the colorbar and manually set the tick labels
col = 'east';
ax = get(axx,'position');
cbar=colorbar(axx,col);
cpos=get(cbar,'Position');
cpos(4)= ax(4)*beta; % Half the thickness
cpos(2)= ax(2); % Move it down
cpos(1) = ax(1)+ax(3)+alpha;
cpos(3) = max(cpos(3)/4,0.0075); % elongate it
set(cbar,'Position',cpos);
set(cbar,'yaxislocation','right')
elseif strcmp(orientation,'custom')
col = 'east';
cbar=colorbar(axx,col);
cpos=get(cbar,'Position');
cpos(1) = cpos(1)+alpha(1) ;
cpos(2) = cpos(2)+alpha(2) ;
cpos(3) = cpos(3)/alpha(3); % elongate it
cpos(4) = cpos(4)/alpha(3); % elongate it
set(cbar,'Position',cpos);
set(cbar,'yaxislocation','right')
else
col = 'east';
cbar=colorbar(axx,col);
cpos=get(cbar,'Position');
cpos(4)=cpos(4)/4; % Half the thickness
%cpos(2)=cpos(2) - 0.11; % Move it down
cpos(1) = cpos(1)+alpha ;
cpos(3) = cpos(3)/3; % elongate it
set(cbar,'Position',cpos,'fontsize',8);
set(cbar,'yaxislocation','right')

end

%set(cbar,'fontsize',6);
if nargin>=4
set(get(cbar,'ylabel'),'string',tit)
end
