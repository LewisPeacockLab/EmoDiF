function [] = compute_horz_indpnt_vol()

total = 0;
for k = 1:3
	total = total + compute_horz_indpnt_vol_per_branch(k); % 6.1676e+16
end

resolution = 0.0001;
x1 = 0+resolution:resolution:1-resolution;
length(x1)^2 * total % 6.1664e+24
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[super_total_entries] = compute_horz_indpnt_vol_per_branch(branch)

resolution = 0.0001;
no_chunks = 8;

switch branch
case 1
	%{
	2.2503e+16
	Branch I: y2 defines the dip and y3 defines the rise
	-1 <= y2 < 0, y2 is the dip so it must fall below zero
	0 < y3 <= 1, y3 is the rise so it must fall above zero
	-1 <= y4 <= y3, y4 can hold any value that is below the rise (y3)
	y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
	%}

	y2 = -1:resolution:1;
	y2 = y2(find(y2 >= -1 & y2 < 0));
	y3 = -1:resolution:1;
	y3 = y3(find(y3 > 0 & y3 <= 1));
	y2_and_y3 = allcomb(y2, y3);
	clear y2; clear y3;

	chunk_idx = [1:size(y2_and_y3,1)/no_chunks:size(y2_and_y3,1);size(y2_and_y3,1)/no_chunks:size(y2_and_y3,1)/no_chunks:size(y2_and_y3,1)];
	for i = 1:size(chunk_idx,2);
		y2_and_y3_chunk = y2_and_y3(chunk_idx(1, i):chunk_idx(2, i), :);
		disp(sprintf('Saving y2_and_y3_chunk_%d', i))
		save(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y2_and_y3_chunk_%d', i), 'y2_and_y3_chunk');
	end
	clear y2_and_y3; clear y2_and_y3_chunk;

	super_total_entries = 0;
	matlabpool close force local; % Closing any old jobs
	matlabpool open;
	for i = 1:size(chunk_idx,2)
		disp(sprintf('Loading y2_and_y3_chunk_%d', i))
		load(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y2_and_y3_chunk_%d.mat', i)); % loads a y2_and_y3_chunk
		total_entries = zeros(size(y2_and_y3_chunk,1), 1);
		parfor j = 1:(size(y2_and_y3_chunk,1))
			y4 = length(-1:resolution:y2_and_y3_chunk(j, 2));
			y1 = length(y2_and_y3_chunk(j, 1)+resolution:resolution:1);
			total_entries(j, 1) = y4 * y1;
		end
		super_total_entries = super_total_entries + sum(total_entries);
       
	end
	matlabpool close;

case 2
	%{
	1.5837e+16
	Branch II: y2 defines the dip and y4 defines the rise
	-1 <= y2 < 0, y2 is the dip so it must fall below zero
	0 < y4 <= 1, y4 is the rise so it must fall above zero
	y2 <= y3 <= y4, y3 can hold any value between the dip and the rise
	y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
	%}

	y2 = -1:resolution:1;
	y2 = y2(find(y2 >= -1 & y2 < 0));
	y4 = -1:resolution:1;
	y4 = y4(find(y4 > 0 & y4 <= 1));
	y2_and_y4 = allcomb(y2, y4);
	clear y2; clear y4;

	chunk_idx = [1:(size(y2_and_y4,1))/no_chunks:(size(y2_and_y4,1));(size(y2_and_y4,1))/no_chunks:(size(y2_and_y4,1))/no_chunks:(size(y2_and_y4,1))];
	for i = 1:size(chunk_idx,2)
		y2_and_y4_chunk = y2_and_y4(chunk_idx(1, i):chunk_idx(2, i), :);
		disp(sprintf('Saving y2_and_y4_chunk_%d', i))
		save(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y2_and_y4_chunk_%d', i), 'y2_and_y4_chunk');
	end
	clear y2_and_y4; clear y2_and_y4_chunk;

	super_total_entries = 0;
	matlabpool close force local; % Closing any old jobs
	matlabpool open;
	for i = 1:size(chunk_idx)
		disp(sprintf('Loading y2_and_y4_chunk_%d', i))
		load(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y2_and_y4_chunk_%d.mat', i)); % loads a y2_and_y4_chunk
		total_entries = zeros((size(y2_and_y4_chunk,1)),1);
		parfor j = 1:(size(y2_and_y4_chunk,1))
			y3 = length(y2_and_y4_chunk(j, 1):resolution:y2_and_y4_chunk(j, 2));
			y1 = length(y2_and_y4_chunk(j, 1)+resolution:resolution:1);
			total_entries(j, 1) = y3 * y1;
		end
		super_total_entries = super_total_entries + sum(total_entries);
	end
	matlabpool close;

case 3
	%{
	2.3336e+16
	Branch III: y3 defines the dip and y4 defines the rise
	-1 <= y3 < 0, y3 is the dip so it must fall below zero
	0 < y4 <= 1, y4 is the rise so it must fall above zero
	y3 < y1 <= 1, y1 can hold any value that is above the dip (y3)
	y3 <= y2 <= 1, y2 can hold any value that is above the dip (y3)
	%}

	y3 = -1:resolution:1;
	y3 = y3(find(y3 >= -1 & y3 < 0));
	y4 = -1:resolution:1;
	y4 = y4(find(y4 > 0 & y4 <= 1));
	y3_and_y4 = allcomb(y3, y4);
	clear y3; clear y4;

	chunk_idx = [1:(size(y3_and_y4,1))/no_chunks:size(y3_and_y4,1);size(y3_and_y4,1)/no_chunks:size(y3_and_y4,1)/no_chunks:size(y3_and_y4,1)];
	for i = 1:size(chunk_idx,2)
		y3_and_y4_chunk = y3_and_y4(chunk_idx(1, i):chunk_idx(2, i), :);
		disp(sprintf('Saving y3_and_y4_chunk_%d', i))
		save(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y3_and_y4_chunk_%d', i), 'y3_and_y4_chunk');
	end
	clear y3_and_y4; clear y3_and_y4_chunk;

	super_total_entries = 0;
	matlabpool close force local; % Closing any old jobs
	matlabpool open;
	for i = 1:size(chunk_idx,2)
		disp(sprintf('Loading y3_and_y4_chunk_%d', i))
		load(sprintf('/Users/tw24955/imdif/analysis/P-CIT/horz_indpt_chance_pp/y3_and_y4_chunk_%d.mat', i)); % loads a y3_and_y4_chunk
		total_entries = zeros(size(y3_and_y4_chunk,1), 1);
		parfor j = 1:size(y3_and_y4_chunk,1)
			y1 = length(y3_and_y4_chunk(j, 1)+resolution:resolution:1);
			y2 = length(y3_and_y4_chunk(j, 1):resolution:1);
			total_entries(j, 1) = y1 * y2;
		end
		super_total_entries = super_total_entries + sum(total_entries);
	end
	matlabpool close;

otherwise, error('Not a valid branch');
end
end

% super_total_entries

function A = allcomb(varargin)

% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
%    of the number of elements of the N inputs. This functionality is also
%    known as the Cartesian Product. The arguments can be numerical and/or
%    characters, or they can be cell arrays.
%
%    Examples:
%       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
%       % -> [ 1  -3   0
%       %      1  -3   1
%       %      1   8   0
%       %        ...
%       %      5  -3   1
%       %      5   8   1 ] ; % a 12-by-3 array
%
%       allcomb('abc','XY') % character arrays
%       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
%
%       allcomb('xy',[65 66]) % a combination
%       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
%
%       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
%       % -> {  'hello'  'Joe'        [99999]
%       %       'hello'  'Joe'             []
%       %       'hello'  [1x3 double] [99999]
%       %       'hello'  [1x3 double]      []
%       %       'Bye'    'Joe'        [99999]
%       %       'Bye'    'Joe'             []
%       %       'Bye'    [1x3 double] [99999]
%       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
%
%    ALLCOMB(..., 'matlab') causes the first column to change fastest which
%    is consistent with matlab indexing. Example: 
%      allcomb(1:2,3:4,5:6,'matlab') 
%      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
%
%    If one of the arguments is empty, ALLCOMB returns a 0-by-N empty array.
%    
%    See also NCHOOSEK, PERMS, NDGRID
%         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)

% Tested in Matlab R2015a
% version 4.1 (feb 2016)
% (c) Jos van der Geest
% email: samelinoa@gmail.com

% History
% 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
%     added option to let the first input run fastest (suggestion by JD)
% 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
%     output by NDGRID. Thanks to Jan Simon, for showing this little trick
% 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
% return an empty output.
% 2.1 (feb 2011). A cell as input argument caused the check on the last
%      argument (specifying the order) to crash.
% 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
% 3.0 (may 2012) removed check for doubles so character arrays are accepted
% 4.0 (feb 2014) added support for cell arrays
% 4.1 (feb 2016) fixed error for cell array input with last argument being
%     'matlab'. Thanks to Richard for pointing this out.

narginchk(1,Inf) ;

NC = nargin ;

% check if we should flip the order
if ischar(varargin{end}) && (strcmpi(varargin{end},'matlab') || strcmpi(varargin{end},'john')),
    % based on a suggestion by JD on the FEX
    NC = NC-1 ;
    ii = 1:NC ; % now first argument will change fastest
else
    % default: enter arguments backwards, so last one (AN) is changing fastest
    ii = NC:-1:1 ;
end

args = varargin(1:NC) ;
% check for empty inputs
if any(cellfun('isempty',args)),
    warning('ALLCOMB:EmptyInput','One of more empty inputs result in an empty output.') ;
    A = zeros(0,NC) ;
elseif NC > 1
    isCellInput = cellfun(@iscell,args) ;
    if any(isCellInput)
        if ~all(isCellInput)
            error('ALLCOMB:InvalidCellInput', ...
                'For cell input, all arguments should be cell arrays.') ;
        end
        % for cell input, we use to indices to get all combinations
        ix = cellfun(@(c) 1:numel(c), args,'un',0) ;
        
        % flip using ii if last column is changing fastest
        [ix{ii}] = ndgrid(ix{ii}) ;
        
        A = cell(numel(ix{1}),NC) ; % pre-allocate the output
        for k=1:NC,
            % combine
            A(:,k) = reshape(args{k}(ix{k}),[],1) ;
        end
    else
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}),[],NC) ;
    end
elseif NC==1,
    A = args{1}(:) ; % nothing to combine

else % NC==0, there was only the 'matlab' flag argument
    A = zeros(0,0) ; % nothing
end
end


