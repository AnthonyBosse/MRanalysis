function TTT = naninterp2(TT,meth)

if nargin<2, meth = 'linear'; end

for k=1:length(TT(1,:))
    if length(find(~isnan(TT(:,k))))>2
       	TTT(:,k) = naninterp(TT(:,k),meth);
        %TTT(:,k) = NaN_interp(TT(:,k),meth); % Ilker, used mine, to avoid linear extrap-- just extend
    else
        TTT(:,k) = TT(:,k);
    end
end
