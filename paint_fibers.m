function out = paint_fibers(args)
% PAINT_FIBERS  Paint stiff vs. floppy fiber regions from coordinates.
%
% Author: Sarah Eldeen (refactor by contributor)
% Date:   2025-09-21
%
% PURPOSE
%   Build masks centered on fiber midpoints and export:
%     - stiff-only masked image
%     - floppy-only masked image
%     - color overlay (stiff=red, floppy=cyan) over grayscale image
%
% INPUT (struct 'args')
%   args.imageFile      : path to the grayscale TIFF (required)
%   args.coordsXlsx     : path to Excel file with columns [x1 y1 x2 y2] (required)
%   args.matPath        : path to MAT file containing logical vectors:
%                         'stiff_indices' and 'floppy_indices' (required)
%   args.indexRange     : [start end] inclusive indices into the rows of coords & indices
%                         (default: [] = use all rows)
%   args.squareSize     : square side length in pixels (default: 50)
%   args.shape          : 'square' (default) or 'disk'
%   args.diskRadius     : if shape='disk', radius in pixels (default: squareSize/2)
%   args.outDir         : output directory (default: same as image)
%   args.savePrefix     : prefix for outputs (default: basename(imageFile))
%   args.verbose        : true/false (default: true)
%
% OUTPUT (struct)
%   .stiffMask, .floppyMask   : logical masks (image size)
%   .stiffImage, .floppyImage : masked images (same class as input)
%   .overlayRGB               : uint8 RGB overlay
%   .usedIdx                  : linear indices (rows) used from inputs
%   .params                   : copy of args (with defaults filled)
%
% REQUIREMENTS
%   MATLAB R2019b+ (for readtable/imread/imwrite)
%
% NOTES
%   - Coordinates are interpreted in pixel units (image column=x, row=y).
%   - Midpoints are computed as [(x1+x2)/2, (y1+y2)/2] per row.
%   - If stiff & floppy overlap, both colors are blended in overlay.

%% -------- Defaults & validation --------
if nargin < 1, args = struct; end
args = withDefaults(args, struct( ...
    'imageFile',   '', ...
    'coordsXlsx',  '', ...
    'matPath',     '', ...
    'indexRange',  [], ...
    'squareSize',  50, ...
    'shape',       'square', ...
    'diskRadius',  [], ...
    'outDir',      '', ...
    'savePrefix',  '', ...
    'verbose',     true ...
));

assert(~isempty(args.imageFile)  && isfile(args.imageFile),  'imageFile not found.');
assert(~isempty(args.coordsXlsx) && isfile(args.coordsXlsx), 'coordsXlsx not found.');
assert(~isempty(args.matPath)    && isfile(args.matPath),    'matPath not found.');
assert(args.squareSize > 0, 'squareSize must be > 0.');

if isempty(args.outDir),  args.outDir = fileparts(args.imageFile); end
if isempty(args.savePrefix)
    [~, base, ~] = fileparts(args.imageFile);
    args.savePrefix = base;
end
if isempty(args.diskRadius), args.diskRadius = args.squareSize/2; end

args.shape = validatestring(args.shape, {'square','disk'});

if ~isfolder(args.outDir), mkdir(args.outDir); end

%% -------- Load inputs --------
if args.verbose, fprintf('Reading image: %s\n', args.imageFile); end
I = imread(args.imageFile);
imgClass = class(I);
I = double(I);                                % work in double internally

if args.verbose, fprintf('Reading coordinates: %s\n', args.coordsXlsx); end
tab = readtable(args.coordsXlsx);
% Expect columns: [x1 y1 x2 y2] starting in col 2 like original script.
% Robustly find first 4 numeric columns.
numTab = tab(:, varfun(@isnumeric, tab, 'OutputFormat','uniform'));
assert(width(numTab) >= 4, 'Expected at least 4 numeric columns (x1 y1 x2 y2).');
numArr = table2array(numTab(:,1:4));
x1 = numArr(:,1); y1 = numArr(:,2);
x2 = numArr(:,3); y2 = numArr(:,4);
xmid = (x1 + x2) / 2;
ymid = (y1 + y2) / 2;
midpoints = [xmid, ymid];

if args.verbose, fprintf('Loading indices: %s\n', args.matPath); end
S = load(args.matPath);
assert(isfield(S,'stiff_indices') && isfield(S,'floppy_indices'), ...
    'MAT must contain logical vectors ''stiff_indices'' and ''floppy_indices''.');

stiff_indices  = S.stiff_indices;
floppy_indices = S.floppy_indices;

% Make sure they are column logicals with same length as midpoints
stiff_indices  = logical(stiff_indices(:));
floppy_indices = logical(floppy_indices(:));
nRows = size(midpoints,1);
assert(numel(stiff_indices)  == nRows, 'stiff_indices length mismatch.');
assert(numel(floppy_indices) == nRows, 'floppy_indices length mismatch.');

% Apply optional index range
if ~isempty(args.indexRange)
    validateattributes(args.indexRange, {'numeric'}, {'numel',2,'integer','positive','increasing'});
    a = max(1, args.indexRange(1));
    b = min(nRows, args.indexRange(2));
    usedIdx = (a:b).';
else
    usedIdx = (1:nRows).';
end

coords  = midpoints(usedIdx, :);
stiff_v = stiff_indices(usedIdx);
floppy_v= floppy_indices(usedIdx);

%% -------- Build masks --------
[h, w] = size(I);
stiffMask  = false(h, w);
floppyMask = false(h, w);

switch args.shape
    case 'square'
        half = args.squareSize / 2;
        for i = 1:size(coords,1)
            xc = coords(i,1); yc = coords(i,2);
            xL = max(1, floor(xc - half));
            xR = min(w, floor(xc - half + args.squareSize - 1));
            yT = max(1, floor(yc - half));
            yB = min(h, floor(yc - half + args.squareSize - 1));
            if stiff_v(i),  stiffMask(yT:yB, xL:xR)  = true; end
            if floppy_v(i), floppyMask(yT:yB, xL:xR) = true; end
        end
    case 'disk'
        R = args.diskRadius;
        [XX, YY] = meshgrid(1:w, 1:h);
        for i = 1:size(coords,1)
            xc = coords(i,1); yc = coords(i,2);
            disk = (XX - xc).^2 + (YY - yc).^2 <= R^2;
            if stiff_v(i),  stiffMask  = stiffMask  | disk; end
            if floppy_v(i), floppyMask = floppyMask | disk; end
        end
end

% Generate masked images
stiffImg  = I;  stiffImg(~stiffMask)   = 0;
floppyImg = I;  floppyImg(~floppyMask) = 0;

% Cast back to original class for saving
stiffImgOut  = cast(imclip(stiffImg, imgClass),  imgClass);
floppyImgOut = cast(imclip(floppyImg, imgClass), imgClass);

%% -------- Build overlay (stiff=red, floppy=cyan) --------
gray = I;
gray = gray - min(gray(:));
if max(gray(:)) > 0, gray = gray ./ max(gray(:)); end  % scale to [0,1]
R = gray; G = gray; B = gray;

% color the masks (alpha-blend)
alpha = 0.6;
R(stiffMask)  = (1-alpha)*R(stiffMask)  + alpha*1.0; % red
G(stiffMask)  = (1-alpha)*G(stiffMask)  + alpha*0.0;
B(stiffMask)  = (1-alpha)*B(stiffMask)  + alpha*0.0;

R(floppyMask) = (1-alpha)*R(floppyMask) + alpha*0.0; % cyan (G+B)
G(floppyMask) = (1-alpha)*G(floppyMask) + alpha*1.0;
B(floppyMask) = (1-alpha)*B(floppyMask) + alpha*1.0;

overlayRGB = uint8(255 * cat(3, R, G, B));

%% -------- Save outputs --------
stiffPath  = fullfile(args.outDir, sprintf('%s_stiff_masked.tif',  args.savePrefix));
floppyPath = fullfile(args.outDir, sprintf('%s_floppy_masked.tif', args.savePrefix));
overlayPath= fullfile(args.outDir, sprintf('%s_overlay.png',       args.savePrefix));

imwrite(stiffImgOut,  stiffPath);
imwrite(floppyImgOut, floppyPath);
imwrite(overlayRGB,   overlayPath);

if args.verbose
    fprintf('Saved: %s\nSaved: %s\nSaved: %s\n', stiffPath, floppyPath, overlayPath);
end

%% -------- Package outputs --------
out = struct();
out.stiffMask   = stiffMask;
out.floppyMask  = floppyMask;
out.stiffImage  = stiffImgOut;
out.floppyImage = floppyImgOut;
out.overlayRGB  = overlayRGB;
out.usedIdx     = usedIdx;
out.params      = args;

end % function

%% ================= Helpers =================
function s = withDefaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
end

function A = imclip(A, targetClass)
% Clip/scale A to the range of targetClass without changing contrast shape.
switch targetClass
    case {'uint8','uint16','uint32'}
        info = intminmax(targetClass);
        A = max(min(A, info.max), info.min);
    case {'int16','int32'}
        info = intminmax(targetClass);
        A = round(max(min(A, info.max), info.min));
    case 'logical'
        A = A > 0;
    otherwise
        % float types: normalize to [0,1] if out of range
        if any(A(:) < 0) || any(A(:) > 1)
            mn = min(A(:)); mx = max(A(:));
            if mx > mn, A = (A - mn) / (mx - mn); end
        end
end
end

function mm = intminmax(cls)
mm.min = double(intmin(cls));
mm.max = double(intmax(cls));
end
