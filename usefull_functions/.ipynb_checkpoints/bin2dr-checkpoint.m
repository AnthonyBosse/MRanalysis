function bin=bin2dr(x,y,d,X,Y,asc,err)
% BIN2D         2D binning of data vectors (MEAN or MVUE)
%
% bin = bin2d(x,y,d,X,Y,asc,err)
% 
% x,y = vectors or matrices of datapoint-positions, corresponding in both
%       size and data order (important)
% d   = vector or matrix of data to be binned (corresponding to x and y)
% X,Y = bin-specifications for the two dimensions separately
%       (default: 10 bin structure over the range of positions).
%       See BUILDGRID on how to specify!
% asc = Logical 1 or 0, on or off, for ASCII plotting of number of
%       datapoints (default = 1)
% err = vector of measurement errors/weights. If given, the MVUE method
%       is used.   
%
% bin = the results in a structural of matrices and vectors:
%       .n              number of points in bins (basis for mean-values)
%       .mean           estimated mean-values
%       .var            error-variance for the mean-values (sigma^2/N-1)
%       .min  .max      min and max values
%       .x,y            mid-points/positions of bins    (vectors)
%       .xg,yg          limits between bins             (vectors)
%
% The matrices are designed for (.x,.y,.mean)-plotting, so the x-dimension
% runs along the matrices' 2nd dimension (along the rows) and the
% y-dimension along the matrices' 1st dimension (down the columns). 
%
% Simple ASCII plotting of the number of datapoints in the bins is done
% to give an overview of number of data in bins:
% '0-9'=0-9 ; '.'=10-50 ; ':'=50-100 ; '*'=100-1000 ; '#'=1000+
%
% See also BUILDGRID BIN3D BIN1D

%Time-stamp:<Last updated on 06/12/10 at 23:05:08 by even@nersc.no>
%Fil:<d:/home/matlab/bin2d.m>       

error(nargchk(3,6,nargin));
NN=10; % default number of bins
if nargin<7 | isempty(err), err=[]; end
if nargin<6 | isempty(asc), asc=logical(1); end
if nargin<5 | isempty(Y),   Y=NN;   end
if nargin<4 | isempty(X),   X=NN;   end

% The full 2D-structure (both XG and YG) needs to be buildt here, based
% on the spans of all data (in case a non-specific bin-spec is given). 
% The XG-bin-specs are then to be sendt down through to BIN1D.
[X,XG]=buildgrid(X,x);          % buildgrid also filters for matrix input
[Y,YG]=buildgrid(Y,y);          % of bin-specs and 

M=length(Y); N=length(X);

%p=nans(size(d)); %(maybe unneccesary spec)      % data-tracker

% the binning for each y-bin
for i=1:M % loop through M rows/along column ( y = i = | )
  if asc, fprintf('Processing row number %2d of %i:  ',i,M); end
  %if i==M,     index=find(YG(i)<=y & y<=YG(i+1) & ~isnan(d));
  %else         index=find(YG(i)<=y & y< YG(i+1)  & ~isnan(d)); end
  if i==M,      index=find(YG(i)<=y & y<=YG(i+1));
  else          index=find(YG(i)<=y & y< YG(i+1));      end
  if isempty(err)
    %[biny,py]=bin1d(x(index),d(index),XG*j);
    biny=bin1dr(x(index),d(index),XG*j,asc);
  else
    %[biny,py]=bin1d(x(index),d(index),XG*j,err(index));
    biny=bin1dr(x(index),d(index),XG*j,asc,err(index))
  end
  bin.n(i,:)     = biny.n;              % number of data in bin
  bin.mean(i,:)  = biny.mean;           % mean value
  bin.median(i,:)  = biny.median;           % median value
  bin.var(i,:)   = biny.var;            % variance
  bin.min(i,:)   = biny.min;            % min
  bin.max(i,:)   = biny.max;            % max
  %p(index)       = (i-1)*N + py;        % data-tracker (bin-number)
end
if asc, fprintf('\n'); end
 
bin.x=X;        bin.y=Y;                        % grid-mid out
bin.xg=XG;      bin.yg=YG;                      % grid-lims out

% if isempty(p)
%   I=[];
% else
%   [j,i]=ind2sub([N M],p);               % data-tracker (array indices)
%   I=sub2ind([M N],i,j); 
% end

% plot if there are no outarguments
if nargout==0
  disp('Creating plot...');
  fig bin2d 4;clf;
  epcolor(bin.x,bin.y,bin.mean);
  set(gca,'xtick',bin.xg,'ytick',bin.yg);
  %contourf(bin.x,bin.y,bin.mean);
  xlabel x; ylabel y;
  tallibing(bin.x,bin.y,bin.n);
  %ecolorbar(bin.mean);
  scanplot;
  colorbar;
end

% DISABLED for speed:
% p     position in the bin-structure of the input data.
%       A vector corresponding to the datavectors (x,y and
%       d), with numbers showing what bin each datapoint is
%       assigned to.  The numbering of bins is running
%       through each row of x-bins, with y increasing.
% I     Like p, but with matrix indices to .n, .mean and
%       .var instead of bin numbers (See below).  
% This means that the bin numbers in p are _not_ matrix indices. The
% matrix indices are given in I (see above). Furthermore, when plotting
% (.x,.y,.mean) the matrix is "viewed" upside down if .y is increasing
% (not unusual).

