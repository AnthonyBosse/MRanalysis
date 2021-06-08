function answer=isvec(varargin)
% ISVEC         vector-object test and -orientation finder
% 
% answer = isvec(varargin)
%
% answer        = A vector of the following numbers, describing which
%                 input object is 
%                       0) no vector
%                       1) a column vector (along dimension 1)
%                       2) a row vector (along dimension 2)
% 
% See also ISSINGLE ISMATRIX ISARRAY ISVEC

%Time-stamp:<Last updated on 06/04/24 at 11:51:34 by even@nersc.no>
%File:</home/even/matlab/evenmat/isvec.m>

for i=1:length(varargin)
  [M,N]=size(varargin{i});
  find([N M]==1);               % in what dimension is vector _not_ 1 wide
  if length(ans)==1
    %x is a vector 
    answer(i)=ans;              % answer is in which dimension the vector...
  else
    answer(i)=0;
  end
end
%answer=logical(answer); % 2 gets converted to logical 1 in ver 6.5 (R13)

