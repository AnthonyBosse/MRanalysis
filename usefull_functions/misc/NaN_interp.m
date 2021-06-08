function xint=NaN_interp(x,extrapflag)
% xint=NaN_interp(x,extrapflag)
% Find the NaN values in x (not matrix), and replaces
% with the interpolated values output is xint, the same size as x.
%
% extrapflag==1 forces nearest-neighbour extrapolation for start and end
% pieces of record...
if nargin<2
    extrapflag = 0 ;
end

if all(isnan(x))
    xint=x;
    
else
    
    flag=0;
    [m,n]=size(x);
    if n>1; x=x'; flag=1; end
    x1=denan(x); y1=find(~isnan(x));
    y2=[1:length(x)]';
    
    if length(x1)<2
        xint=NaN.*x;
    else
        xint=interp1(y1,x1,y2);
        if extrapflag | strcmp(extrapflag,'extrap')
            id = isnan(xint) ;
            xint=interp1(y2(~id),xint(~id),y2,'nearest','extrap');
        end
    end
    %make sure x and xint are of same dimension
    if flag; xint=xint';end
end
