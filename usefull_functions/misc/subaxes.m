function varargout = subaxes(numRows,numCols,varargin)
%
% subaxes.m--Creates multiple axes laid out in a rectangular grid in the
% current figure.
%
% Somewhat like Matlab's subplot.m, but the spacing between the axes can be
% specified by the user.
%
% The user specifies the numbers of rows and columns of axes. If the
% rectangular grid of numRows rows and numCols columns results in more axes
% than are actually required, specify the number of axes to make with a
% third input argument, numAxes.
%
% The remaining optional input arguments are specified using
% parameter/value pairs. The parameters and their uses are as follows:
% 
%    'RowSpacing' gives the amount of spacing between rows of axes, and is
%    specified in units of axis heights (i.e., a RowSpacing of 0.5 results
%    in half an axis height of spacing between rows of axes). Negative
%    spacing results in overlapping axes; a spacing of NaN causes the
%    default value of the spacing to be used. 
%
%    'ColSpacing' gives the amount of spacing between columns of axes, and is
%    specified in units of axis widths (i.e., a ColSpacing of 0.5 results
%    in half an axis width of spacing between columns of axes). Negative
%    spacing results in overlapping axes; a spacing of NaN causes the
%    default value of the spacing to be used. 
%
%    'PageMargins' allows the user to specify the width of the margins
%    between the sides of the page and the axes. By default, subaxes.m uses
%    margins of [0.13 0.095 0.11 0.075] ([left right bottom top] in
%    normalised units), the same margins used by Matlab's subplot.m.
%    Specify your own set of margins with a different 4-element vector.
%    Alternatively, use the string 'tight' to use the smallest practical
%    margins (depending on your printer, the values for "tight" margins may
%    have to be changed).       
%
% Subaxes.m returns AxHndls, a vector of handles to the axes created, and
% row and col, two vectors giving the row and column number of each set of
% axes.  Optionally, subaxes.m also returns Boolean vectors IsTop,
% IsBottom, IsLeft and IsRight, which indicate whether each subplot is in
% the topmost row, the bottommost row, the leftmost column and/or the
% rightmost column. These vectors are useful for determining which set of
% axes should receive xlabels and titles, etc.
%
% Syntax: [<AxHndls>,<row,col,IsTop,IsBottom,IsLeft,IsRight>] = subaxes(numRows,numCols,<numAxes>,<'RowSpacing',RowSpacing>,<'ColSpacing',ColSpacing>,<'PageMargins',PageMargins>)
%
% e.g.,   [AxHndls,row,col,IsTop,IsBottom,IsLeft,IsRight] = subaxes(3,4)
% e.g.,   [AxHndls,row,col,IsTop,IsBottom,IsLeft,IsRight] = subaxes(3,4,10,'RowSpacing',0,'PageMargins','tight')
% e.g.,   subaxes(3,4,'RowSpacing',0.1,'PageMargins',[0.13 0.095 0.2 0.2])
% e.g.,   subaxes(3,4,11,'RowSpacing',0.1,'ColSpacing',0.1,'PageMargins',[0.13 0.095 0.4 0.4])

% Developed in Matlab 6.5.0.180913a (R13) on Win98.
% Kevin Bartlett(kpb@hawaii.edu), 2003/02/25, 11:40
%------------------------------------------------------------------------------

% N.B., normalised units are used throughout this program.
if ~ismember(nargin,[2:9])
   error([mfilename '.m--Incorrect number of input arguments.']);
end % if

if nargout ~= 0 &  nargout ~= 1 &  nargout ~= 7
   error([mfilename '.m--Incorrect number of output arguments.']);
end % if

% Handle input arguments.
[regArgs,propVals] = parseparams(varargin);

if isempty(regArgs)
   numAxes = numRows * numCols;
else
   numAxes = regArgs{1};
end % if

DefaultVals.RowSpacing = NaN; 
DefaultVals.ColSpacing = NaN; 
DefaultVals.PageMargins = 'default';
allowNewFields = 0;
ParValStruct = parvalpairs(DefaultVals,propVals{:},allowNewFields);
RowSpacing = ParValStruct.RowSpacing;
ColSpacing = ParValStruct.ColSpacing;
PageMargins = ParValStruct.PageMargins;

DefaultMargins = [0.13 0.095 0.11 0.075];
TightMargins = [0.065 0.047 0.055 0.0375];

if strcmp(get(gcf,'PaperOrientation'),'landscape')
   TightMargins = [TightMargins(3) TightMargins(4) TightMargins(1) TightMargins(2) ];
end % if

if isstr(PageMargins)
   
   if strmatch(lower(PageMargins),'default')
      margins = DefaultMargins;
   elseif strmatch(lower(PageMargins),'tight')
      margins = TightMargins;
   else
      error([mfilename '.m--Unrecognised input string.']);
   end % if
   
else
   margins = PageMargins;
end %if

Lmargin = margins(1);
Rmargin = margins(2);
Bmargin = margins(3);
Tmargin = margins(4);

% If no value specified for either or both spacings, use defaults. Defaults
% are chosen so that subaxes.m will produce axes of the similar dimensions as
% produced by Matlab's subplot.m.
if isnan(ColSpacing)
   if numCols>1
      ColSpacing = 0.3210;
   else
      ColSpacing = 0;
   end % if
end % if

if isnan(RowSpacing)
   if numRows>1
      RowSpacing = 0.3208;
   else
      RowSpacing = 0;
   end % if
end % if

% Check for good input values.
if numRows < 1 | numCols < 1 | round(numRows)~=numRows | round(numCols)~=numCols  
   error([mfilename '.m--The numbers of columns and rows must be positive integers.']);
end % if

if RowSpacing < -1 | ColSpacing < -1
   error([mfilename '.m--Axes cannot overlap by greater than 100%. RowSpacing and ColumnSpacing must be at least -1.']);
end % if
   
if numAxes > numRows*numCols | round(numAxes)~=numAxes 
   error([mfilename '.m--The number of axes requested must be an integer no greater than the number of rows times the number of columns.']);
end % if

if numAxes < 1
   error([mfilename '.m--Cannot request less than 1 set of axes.']);
end % if

% Get the figure dimensions.
FigHndl = gcf;
OrigUnits = get(FigHndl,'units');
set(FigHndl,'units','normalized');
FigPos = get(FigHndl,'position');

% Determine axes dimensions.
AxWidth = (1 - (Lmargin+Rmargin)) / (numCols + numCols*ColSpacing - ColSpacing);
AxHeight = (1 - (Bmargin+Tmargin)) / (numRows + numRows*RowSpacing - RowSpacing);

% Convert row and column spacing from fractional to normalised (relative to
% figure) units.
ColSpacingNorm = ColSpacing * AxWidth;
RowSpacingNorm = RowSpacing * AxHeight;

% Determine axes' positions.
ColPos = Lmargin + [0:(numCols-1)].*(AxWidth+ColSpacingNorm);
RowPos = Bmargin + [0:(numRows-1)].*(AxHeight+RowSpacingNorm);

% Want to number axes from top to bottom.
RowPos = fliplr(RowPos);

% Create a matrix describing the axes layout, with ones showing where a set
% of axes is to exist, zeroes elsewhere.
axMat = zeros(numRows,numCols);
axMat = axMat';
axMat(1:numAxes) = 1;
axMat = axMat';

% Create axes.
AxHndls = NaN * ones(1,(numRows*numCols));
row = AxHndls;
col = AxHndls;
IsTop = AxHndls;
IsBottom = AxHndls;
IsLeft = AxHndls;
IsRight = AxHndls;

AxCount = 1;

for RowCount = 1:numRows
   
   for ColCount = 1:numCols

      if AxCount > numAxes         
         break;
      end % if
    
      AxHndls(AxCount) = axes('units','normalized','position',[ColPos(ColCount) RowPos(RowCount) AxWidth AxHeight],'box','on');
      row(AxCount) = RowCount;
      col(AxCount) = ColCount;
      IsTop(AxCount) = (RowCount == 1);
      
      if (RowCount == numRows) | axMat(RowCount+1,ColCount)==0
         IsBottom(AxCount) = 1;
      else
         IsBottom(AxCount) = 0;
      end
   
      IsLeft(AxCount) = (ColCount == 1);

      if (ColCount == numCols) | axMat(RowCount,ColCount+1)==0
         IsRight(AxCount) = 1;
      else
         IsRight(AxCount) = 0;
      end
      
      AxCount = AxCount + 1;
      
   end % for
   
end % for

AxHndls = AxHndls(1:numAxes);
row = row(1:numAxes);
col = col(1:numAxes);
IsTop = IsTop(1:numAxes);
IsBottom = IsBottom(1:numAxes);
IsLeft = IsLeft(1:numAxes);
IsRight = IsRight(1:numAxes);

set(FigHndl,'units',OrigUnits);

% Matlab will sometimes move a title up or down when the y-limits of one of 
% these axes is changed. Work-around: create empty titles.m with units set 
% to 'normalized'.
titleHndls = get(AxHndls,'title');

if iscell(titleHndls)
    titleHndls = cat(1,titleHndls{:});
end % if

set(titleHndls,'string','','units','normalized');

% Handle output arguments.
if nargout == 1
   varargout{1} = AxHndls;
elseif nargout == 7;
   varargout{1} = AxHndls;
   varargout{2} = row;
   varargout{3} = col;
   varargout{4} = IsTop;
   varargout{5} = IsBottom;
   varargout{6} = IsLeft;
   varargout{7} = IsRight;
end % if

