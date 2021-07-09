function y = no_mean_median(x,dim,flag)
%MEDIAN Median value.
%   For vectors, MEDIAN(x) is the median value of the elements in x.
%   For matrices, MEDIAN(X) is a row vector containing the median value
%   of each column.  For N-D arrays, MEDIAN(X) is the median value of the
%   elements along the first non-singleton dimension of X.
%
%   MEDIAN(X, DIM) takes the median along the dimension DIM of X.
%
%   MEDIAN(..., MISSING) specifies how NaN (Not-A-Number) values
%   are treated. The default is 'includenan':
%
%   'includenan' - the median of a vector containing NaN values is also NaN.
%   'omitnan'    - the median of a vector containing NaN values is the
%                  median of all its non-NaN elements. If all elements
%                  are NaN, the result is NaN.
%
%   Example: If X = [1 2 4 4; 3 4 6 6; 5 6 8 8; 5 6 8 8];
%
%   then median(X) is [4 5 7 7] and median(X,2) is [3; 5; 7; 7]
%
%   Class support for input X:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also MEAN, STD, MIN, MAX, VAR, COV, MODE.

%   Copyright 1984-2015 The MathWorks, Inc.

if isempty(x) % Not checking flag in this case
    if nargin == 1 || (nargin == 2 && ischar(dim))
        
        % The output size for [] is a special case when DIM is not given.
        if isequal(x,[])
            if isinteger(x)
                y = zeros('like',x);
            else
                y = nan('like',x);
            end
            return;
        end
        
        % Determine first nonsingleton dimension
        dim = find(size(x)~=1,1);
        
    end
    
    s = size(x);
    if dim <= length(s)
        s(dim) = 1;                  % Set size to 1 along dimension
    end
    if isinteger(x)
        y = zeros(s,'like',x);
    else
        y = nan(s,'like',x);
    end
    
    return;
end

omitnan = false;
dimSet = true;
if nargin == 1
    dimSet = false;
elseif nargin == 2 && ischar(dim)
    flag = dim;
    dimSet = false;
end

if nargin == 2 && dimSet == false || nargin == 3
    len = max(length(flag), 1);
    
    if ~isrow(flag)
        error(message('MATLAB:median:unknownFlag'));
    end
    
    s = strncmpi(flag, {'omitnan', 'includenan'}, len);
    
    if ~any(s)
        error(message('MATLAB:median:unknownFlag'));
    end
    
    omitnan = s(1);
end

sz = size(x);

if dimSet
    if ~isscalar(dim) || ~isreal(dim) || floor(dim) ~= ceil(dim) || dim < 1
        error(message('MATLAB:getdimarg:dimensionMustBePositiveInteger'));
    end
    
    if dim > numel(sz)
        y = x;
        return;
    end
end


if isvector(x) && (~dimSet || sz(dim) > 1)
    % If input is a vector, calculate single value of output.
    x = sort(x);
    nCompare = length(x);
    if isnan(x(nCompare))        % Check last index for NaN
        if omitnan
            nCompare = find(~isnan(x), 1, 'last');
            if isempty(nCompare)
                y = nan('like',x([])); % using x([]) so that y is always real
                return;
            end
        else
            y = nan('like',x([])); % using x([]) so that y is always real
            return;
        end
    end
    half = floor(nCompare/2);
    y = x(half+1);
    if 2*half == nCompare        % Average if even number of elements
        y = max(x(half),y);
    end
else
    if ~dimSet              % Determine first nonsingleton dimension
        dim = find(sz ~= 1,1);
        
    end
    
    % Reshape and permute x into a matrix of size sz(dim) x (numel(x) / sz(dim))
    if dim ~= 1
        n1 = prod(sz(1:dim-1));
        n2 = prod(sz(dim+1:end));
        if n2 == 1 % because reshape(x, [n1, n2, 1]) errors for sparse matrices
            x = reshape(x, [n1, sz(dim)]); 
        else
            x = reshape(x, [n1, sz(dim), n2]);
        end
        x = permute(x, [2 1 3]);
    end

    x = reshape(x,size(x,1),[]);
    % Sort along columns
    x = sort(x, 1);
    if ~omitnan || all(~isnan(x(end, :)))
        % Use vectorized method with column indexing.  Reshape at end to
        % appropriate dimension.
        nCompare = sz(dim);          % Number of elements used to generate a median
        half = floor(nCompare/2);    % Midway point, used for median calculation

        y = x(half+1,:);
        if 2*half == nCompare
            y = max(x(half,:),y);
        end

        y(isnan(x(nCompare,:))) = NaN;   % Check last index for NaN

    else

        % Get median of the non-NaN values in each column.
        y = nan(1, size(x, 2), 'like', x([])); % using x([]) so that y is always real

        % Number of non-NaN values in each column
        n = sum(~isnan(x), 1);

        % Deal with all columns that have an odd number of valid values
        oddCols = find((n>0) & rem(n,2)==1);
        oddIdxs = sub2ind(size(x), (n(oddCols)+1)/2, oddCols);
        y(oddCols) = x(oddIdxs);

        % Deal with all columns that have an even number of valid values
        evenCols = find((n>0) & rem(n,2)==0);
        evenIdxs = sub2ind(size(x), n(evenCols)/2, evenCols);
        y(evenCols) = max( x(evenIdxs), x(evenIdxs+1) );

    end
    % Now reshape output.
    sz(dim) = 1;
    y = reshape(y, sz);
end

