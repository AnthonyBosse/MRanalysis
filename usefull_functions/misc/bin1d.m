function bin=bin1d(x,d,X,asc,err)
% BIN1D         1D binning of data vectors (MEAN or MVUE)
%
% bin = bin1d(x,d,X,asc,err)
% 
% x   = Vector or matrix of datapoint-positions corresponding to ...
% d   = vector or matrix of data to be binned.
% X   = Bin-specification (default: 10-bin-structure over the range of
%       positions). See BUILDGRID for description of input options! 
% asc = Logical 1 or 0, on or off, for ASCII plotting of number of
%       datapoints (default = 1)
% err = vector of measurement errors/weights. If given, the MVUE method
%       is used.   
%
% bin = the results in a structural of vectors:
%       .n              number of points in bin (basis for mean-values)
%       .mean           estimated mean-values 
%       .var            variance (sigma^2) for the bin. To calculate
%                       error-variance for the mean-values
%                       (sigma^2/N-1), use .n as N.  
%       .min  .max      min and max values 
%       .x              mid-points/positions of bins	(vector)
%       .xg             limits between bins		(vector)
%
% Simple ASCII plotting of the number of datapoints in the bins is done
% to give an overview of number of data in bins:
% '0-9'=0-9 ; '.'=10-50 ; ':'=50-100 ; '*'=100-1000 ; '#'=1000+
%
% See also BUILDGRID BIN2D BIN3D BINSUB

%Time-stamp:<Last updated on 09/03/05 at 15:36:25 by even@nersc.no>
%File:</Users/even/matlab/evenmat/bin1d.m>       

%error(nargchk(2,4,nargin));
NN=10; % default number of bins
%if nargin<5 | isempty(opt), opt=''; end
if nargin<5 | isempty(err), err=[]; end
if nargin<4 | isempty(asc), asc=logical(1); end
if nargin<3 | isempty(X),   X=NN;   end

% build X-bin structure
[X,XG]=buildgrid(X,x);		% buildgrid also filters for matrix input
                                % of bin-specs
N=length(X);

%initialization of variables (only here in 1D). BIN2D and -3D
[bin.n,bin.mean,bin.var,bin.min,bin.max,bin.prc1,bin.prc5,bin.prc25,bin.prc50,bin.prc75,bin.prc95,bin.prc99]=deal(nan(1,N));
%p=nans(size(d));                    % data-tracker

% the binning _and_ calculations for each x-bin
for j=1:N % loop through N "columns"/along row ( x = j = - )
  if j==N,	index=find(XG(j)<=x & x<=XG(j+1) & ~isnan(d));  
  else		index=find(XG(j)<=x & x< XG(j+1) & ~isnan(d));	end
  bin.n(j) = length(index);             % Number of data in this bin 
  %p(index) = j;                     % Data-tracker
  if bin.n(j)~=0 
    if isempty(err)	%%%%%%%%%%%% REGULAR AVERAGE %%%%%%%%%%%%%%%%%
  %  xx = robust_mean(d(index),1);
      bin.mean(j)=nanmean(d(index));       % (these are the CORE bin-mean-
      if bin.n(j)==1                    %  calculations for 
        bin.var(j) =0;		%  the BINXD functions)
        bin.min(j) =d(index);
        bin.max(j) =d(index);
        %bin.var(j)=err(index); % using the instrument error
      else      
        bin.var(j)   = var(d(index)); 
	%bin.var(j)   = var(d(index))/(bin.n(j)-1); 
	[bin.min(j) bin.max(j)] = mima(d(index));
      end

	%%%% get some percentiles, Anthony
	bin.prc1(j)=prctile(d(index),[1]);
	bin.prc5(j)=prctile(d(index),[5]);
	bin.prc25(j)=prctile(d(index),[25]);
	bin.prc50(j)=prctile(d(index),[50]);
	bin.prc75(j)=prctile(d(index),[75]);
	bin.prc95(j)=prctile(d(index),[95]);
	bin.prc99(j)=prctile(d(index),[99]);
    else		%%%%%%%%%%%% MVUE AVERAGE %%%%%%%%%%%%%%%%%%%%
keyboard
[bin.mean(j),bin.var(j)]=mvue(d(index),err(index));
      if bin.n(j)==1  
	bin.min(j)=NaN; bin.max(j)=NaN; 
      else
	[bin.min(j) bin.max(j)] = mima(d(index));
      end
    end			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  % print a primitive grid in Command Window
  if asc
    if N<60,      ASCIIplot(bin.n(j)); 
    elseif j==1,  fprintf('[row too long to show]'); end  
  end
end
if asc, fprintf('\n'); end

bin.x  = X;					% grid-mid out
bin.xg = XG;					% grid-lims out

% plot if there are no outarguments
if nargout==0
  disp('Creating plot...');
  figure% bin1d 4;clf;
  %plot(bin.x,bin.mean,'o-');
  errorbar(bin.x,bin.mean,sqrt(bin.var./(bin.n-1)),'o-');
  grid on
  %  tallibing(bin.x,bin.mean,bin.n,[],'k');
end

%----------------------------------------------------------
function [mean,err] = mvue(d,sigma)
% calculate the mvue-mean and -error
div=sum(sigma.^(-2));	% divisor
w=sigma.^(-2)/div;	% weight (MVUE)
mean=d(:)'*w(:);	% mean=sum(d_i*w_i)
err=1/div;		% error-variance
%----------------------------------------------------------
function ASCIIplot(n)
% plot an ASCII symbol for different values 
if      n<10,   fprintf('%1i',n);
elseif  n<50,   fprintf('.');
elseif  n<100,  fprintf(':');
elseif  n<1000, fprintf('*');
else            fprintf('#');
end
%----------------------------------------------------------



% DISABLED for speed:
% p     position in the bin-structure of the input data.
%       A vector corresponding to the datavectors (x,y and
%       d), with numbers showing what bin each datapoint is
%       assigned to.  The numbering of bins is running
%       through the bins with increasing x. For
%       1D-binning, this is the index to the bin-vectors.

