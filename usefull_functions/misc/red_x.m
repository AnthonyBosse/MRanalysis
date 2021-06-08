function B = red_x(A)
if length(size(A)) == 3
B = (A(:,:,1:end-1)+A(:,:,2:end))/2;
else
B = (A(:,1:end-1)+A(:,2:end))/2;
end

