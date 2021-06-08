function [X,XG]=buildgrid(X,x)
% BUILDGRID     builds two vectors describing a 1D-bin structure
%
% [X,XG] = buildgrid(X,x)
%
% Input:
% x  = vector or matrix of the positions of the data that are to be
%      binned
% X  = bin-specifications interpreted as follows:
% 
%      single real number       =  the number of bins (N)
%      single imaginary number  =  the bin-size (dx)
%      real valued vector       =  the bins' mid-point positions
%      imaginary valued vector  =  the limits of the bins
%
%      The first two options result in an evenly spaced grid over the range
%      of x (here x needs to be given).  The last two options gives the
%      bin-structure explicitly either by its mid-point positions or by its
%      limits, and the values of the other is distributed between them. This
%      is useful for making inhomogeneous bin-structures like
% 
%      |.| . |  .  |   .   |    .    |   (dots=positions; lines=limits)
%
%      The imaginary input is simply used as an option-flag. Just multiply
%      your input value by the imaginary unit i.  If X is given as a matrix,
%      the matrix must be invariate in one dimension. Then the other
%      dimension, along which the values change, is used as bin-specs.
%
%      Matrix input is only possible with regular-position-matrices
%      (matrices with repeated rows or columns).
%
% X  = length N vector of the bins' mid-point positions
% XG = length N+1 vector of the limits of the bins
%
% See also BIN1D BIN2D BIN3D MAT2VEC

%Time-stamp:<Last updated on 06/04/24 at 13:19:50 by even@nersc.no>
%File:</home/even/matlab/evenmat/buildgrid.m>

%error(nargchk(1,2,nargin));		% INPUT-TESTS:
if isempty(X)
  X=20;
end
%x=x(:);				% secure that x is a vector
X=mat2vec(X);			% secure that X is a vector

if issingle(X)
  if ~exist('x') | ~(isvec(x) | ismatrix(x) | isarray(x))
    error(['Single valued bin-specification needs to be ',...
	   'followed by a vector or regular-grid-matrix ',...
	   'of data-point positions!']);
  else
    r=mima(x);
    if isreal(X),	N=X;	dx=diff(r)/N;		% dx
    else			dx=imag(X);
    end   
    XG=r(1):dx:r(2);					% XG
    if r(2)>XG(end),	XG=[XG XG(end)+dx];	end
    X=findpoints(XG);					% X
  end
else
  if ismatrix(X),	X=mat2vec(X);		end
  if isreal(X),		XG=findbins(X);
  else			XG=imag(X);		X=findpoints(XG);
  end
end


%---------------------------------------------------------------------
function X=findpoints(XG)
% X = points half way between all values in XG
X=uplus(XG(1:end-1))+uplus(diff(XG)/2); % preserves shape
%---------------------------------------------------------------------
function XG=findbins(X)
% XG = points half way between all values of X _and_ outside the
% endpoints. The outer limits have same distance from X's endpoints as
% the limits just inside.
dim=isvec(X);			% preserves shape
diff(X)/2; cat(dim,ans,ans(end));%[ans ans(end)];
XG=uplus(X)+uplus(ans);
XG=cat(dim,X(1)-(XG(1)-X(1)),XG);%[X(1)-(XG(1)-X(1)) XG];


