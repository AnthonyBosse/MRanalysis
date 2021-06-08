function parValStruct = parvalpairs(varargin)
%
% parvalpairs.m--Fills the structure parValStruct with parameter/value
% pairs. (Based on Firing Group's fillstruct.m).
%
% Parvalpairs.m can be used as a client function for handling input
% arguments to a parent function. For example, imagine that you want the
% function myfunction.m to accept input in the form of parameter/value
% pairs passed to its varargin variable as a comma-separated list:
%
%                    myfunction('a',1,'b',2,'c',3)
%
% Within myfunction.m, the parameters 'a', 'b' and 'c', and their
% respective values can be extracted from the varargin variable with the
% command:
%
%             parValStruct = parvalpairs(varargin{:})
%
% Default values can be defined and over-ridden inside the parent function
% as follows:
%
%    defaultVals.a = -99; defaultVals.b = -99;
%    parValStruct = parvalpairs(defaultVals,varargin{:})
%
% In the example above, a parameter 'c' was returned as a field of the
% parValStruct variable, even though no default value was defined for it.
% Parvalpairs.m will create as many parameter/value pairs as the user cares
% to pass to the varargin variable. To prevent this, use the optional input
% argument ctrlVec, a vector of integers that controls various aspects of
% parvalpair.m's behaviour. Set the first integer of ctrlVec to 0 to
% prevent the creation of new fields, or 1 to allow it.
%
% The second integer in the ctrlVec vector controls case-sensitivity. By
% default, parvalpairs.m treats 'thisFieldName' and 'thisfieldname' as
% different field names. With the second integer of ctrlVec set to 0,
% parvalpairs.m will be case-INsensitive when it comes to fieldnames.
%
% Syntax: parValStruct = parvalpairs(<defaultStruct>,par1,val1,par2,val2,...,<[allowNewFields isCaseSensitive]>)
%
% e.g.,   defaultStruct.a = pi; defaultStruct.b = 'hello'; defaultStruct.c = [1;2;3];
%         parValStruct = parvalpairs(defaultStruct,'a',pi/2,'b','bye',0)
%
% e.g.,   defaultStruct.a = pi; defaultStruct.b = 'hello'; defaultStruct.c = [1;2;3];
%         allowNewFields=1; isCaseSensitive=0; ctrlVec=[allowNewFields isCaseSensitive];
%         parValStruct = parvalpairs(defaultStruct,'a',pi/2,'b','bye','z','new field',ctrlVec)
%
% e.g.,   parValStruct = parvalpairs('Position',[256 308 512 384],'Units','pixels','Color',[1 0 1])

% Developed in Matlab 6.1.0.450 (R12.1) on Linux. Kevin
% Bartlett(kpb@hawaii.edu), 2003/04/08, 11:49
%--------------------------------------------------------------------------

% Handle input arguments.
if nargin == 0
    %error([mfilename '.m--Incorrect number of input arguments.']);
    parValStruct = struct([]);
    return;
end % if

% Extract existing structured variable from input arguments, if it exists.
%defaultStruct = [];

% ...2005-04-08: Use empty structure, rather than [], in order to permit
% use of "fieldnames" command later.
defaultStruct = struct([]);

if isstruct(varargin{1})

    defaultStruct = varargin{1};

    if nargin > 1
        varargin = varargin(2:end);
    else
        varargin = {};
    end % if

end % if

parValStruct = defaultStruct;

% If a control vector has been specified, there will be an odd number
% of input arguments remaining.
isCaseSensitive_default = 1;
allowNewFields_default = 1;

if rem(length(varargin),2) == 1

    ctrlVec = varargin{end};
    varargin = varargin(1:end-1);

    % Allow for backwards compatiblity by permitting ctrlVec to be of length 1,
    % containing only the allowNewFields variable.
    if length(ctrlVec) == 1
        isCaseSensitive_default = 1;
        ctrlVec = [ctrlVec isCaseSensitive_default];
    end % if

    allowNewFields = ctrlVec(1);
    isCaseSensitive = ctrlVec(2);

else
    allowNewFields = allowNewFields_default;
    isCaseSensitive = isCaseSensitive_default;
end % if

if ~ismember(allowNewFields,[1 0])
    error([mfilename '.m--Value for "allowNewFields" must be 1 or 0.']);
end % if

if ~ismember(isCaseSensitive,[1 0])
    error([mfilename '.m--Value for "isCaseSensitive" must be 1 or 0.']);
end % if

% Remaining input arguments should be parameter/value pairs.
existingFieldNames = fieldnames(defaultStruct);
lowerExistingFieldNames = lower(existingFieldNames);

for ArgCount = 1:2:length(varargin)

    CurrFieldName = varargin{ArgCount};

    if ~isstr(CurrFieldName)
        error([mfilename '.m--Parameter names must be strings.']);
    end % if

    CurrField = varargin{ArgCount+1};

    % Find out if field already exists.
    if isCaseSensitive == 1

        if ismember(CurrFieldName,existingFieldNames)
            fieldExists = 1;
            fieldNameToInsert = CurrFieldName;
        else
            fieldExists = 0;
            fieldNameToInsert = CurrFieldName;
        end % if

    else

        if ismember(lower(CurrFieldName),lowerExistingFieldNames)
            fieldExists = 1;
            matchIndex = strmatch(lower(CurrFieldName),lowerExistingFieldNames,'exact');
            fieldNameToInsert = existingFieldNames{matchIndex};
        else
            fieldExists = 0;
            fieldNameToInsert = CurrFieldName;
        end % if

    end % if

    % If new fields are not permitted to be added, test that field already exists.
    if allowNewFields == 0 & fieldExists == 0
        
        % If the default structure is empty, and the user is not permitting
        % the addition of new fields, there is no point in running this
        % program; probably the user doesn't intend this.
        if isempty(defaultStruct)
            error([mfilename '.m--Need to permit the addition of new fields if no default structure specified']);
        end % if

        if isCaseSensitive == 1
            error([mfilename '.m--Attempt to add new parameter to existing set (case-sensitive). Use allowNewFields=1 to allow this.']);
        else
            error([mfilename '.m--Attempt to add new parameter to existing set. Use allowNewFields=1 to allow this.']);
        end % if

    end % if

    % 2005-04-08--setfield is now deprecated. In addition, now using empty
    % structure instead of [] as starting point, which causes an error.
    % Work around this.
    %parValStruct = setfield(parValStruct,fieldNameToInsert,CurrField);

    if isempty(parValStruct)
        parValStruct = struct(fieldNameToInsert,CurrField);
    else
        parValStruct.(fieldNameToInsert) = CurrField;
    end % if
            
end % for

