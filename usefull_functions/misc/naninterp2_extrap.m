function TTT = naninterp2_extrap(TT)
for k=1:length(TT(1,:))
    if length(find(~isnan(TT(:,k))))>2
       	%TTT(:,k) = naninterp(TT(:,k));
        TTT(:,k) = NaN_interp(TT(:,k),1); % Ilker, used mine, to avoid linear extrap-- just extend
    else
        TTT(:,k) = TT(:,k);
    end
end
