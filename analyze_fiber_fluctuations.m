function out = analyze_fiber_fluctuations(args)
% ANALYZE_FIBER_FLUCTUATIONS
% Calculate minimum frequency / positional fluctuations from Olympus OIR line
% scans by tracking fiber centroids via Gaussian fits, then summarizing the
% distribution width (std dev in microns).
%
% Project: Fiber fluctuation
% Author: Sarah Eldeen (refactor by contributor)
% Date: 2025-10-29
%
% USAGE
%   out = analyze_fiber_fluctuations(struct( ...
%       'dataDir',        'path/to/line/scans', ...
%       'pix2um',         1/38.6175, ...     % microns per pixel
%       'downsample',     10, ...            % take every Nth line
%       'roiHalfWidthPx', 13, ...            % half-width around peak (in px)
%       'gofThreshold',   0.80, ...          % R^2 cutoff for valid fit
%       'showExample',    true, ...          % quick preview figure
%       'savePrefix',     '' ...             % e.g. 'fiber_run1' to save CSV/MAT
%   ));
%
% REQUIREMENTS
%   - MATLAB R2020b+ recommended
%   - Curve Fitting Toolbox (for 'fit' with 'gauss1')
%   - Bio-Formats for MATLAB (bfopen): 
%     https://www.mathworks.com/matlabcentral/fileexchange/129249-bioformats-image-toolbox
%
% INPUT (args struct)
%   dataDir        : folder containing .oir files
%   pix2um         : microns per pixel (scalar)
%   downsample     : (int) keep every Nth scanline
%   roiHalfWidthPx : (int) half window size around peak for fitting
%   gofThreshold   : (0..1) R^2 cutoff for accepting a Gaussian fit
%   showExample    : (bool) show a preview image + axes with physical units
%   savePrefix     : (char) if nonempty, saves results as CSV + MAT
%
% OUTPUT (out struct)
%   .files               : file names
%   .linePeriodSeconds   : 1xN seconds per line for each file (best-effort from metadata; NaN if not found)
%   .centroid_um         : 1xN cell of centroid time-series (microns) after GOF filtering
%   .gof                 : 1xN cell of R^2 values (NaN where rejected)
%   .stdDevs_um          : 1xN vector of Gaussian-fitted histogram std devs (microns)
%   .histFits            : 1xN cell of fit objects (or [])
%   .params              : copy of args
%
% NOTES
%   - Metadata key for line speed varies across acquisitions. This script tries
%     a few likely keys and falls back to NaN. If NaN, time axis in the preview
%     figure shows "lines" instead of seconds.
%   - The previous version appended the same plane repeatedly; this fixes that
%     by iterating planes correctly.
%   - The previous version mixed up variables in the example plot; this fixes it.

%% --------- Defaults & validation ---------
if nargin < 1 || isempty(args), args = struct; end
args = withDefaults(args, struct( ...
    'dataDir',        pwd, ...
    'pix2um',         1/38.6175, ...
    'downsample',     10, ...
    'roiHalfWidthPx', 13, ...
    'gofThreshold',   0.80, ...
    'showExample',    true, ...
    'savePrefix',     '' ...
));

mustExistFolder(args.dataDir);
assert(isscalar(args.pix2um) && args.pix2um > 0, 'pix2um must be a positive scalar.');
assert(args.downsample >= 1 && mod(args.downsample,1)==0, 'downsample must be a positive integer.');
assert(args.roiHalfWidthPx >= 1 && mod(args.roiHalfWidthPx,1)==0, 'roiHalfWidthPx must be a positive integer.');
assert(args.gofThreshold >= 0 && args.gofThreshold <= 1, 'gofThreshold must be in [0,1].');

% Toolboxes
assert(exist('bfopen','file')==2, 'Bio-Formats "bfopen" not found on path.');
assert(license('test','Curve_Fitting_Toolbox')==1, 'Curve Fitting Toolbox is required.');

%% --------- Discover files ---------
fileList = dir(fullfile(args.dataDir, '*.oir'));
if isempty(fileList)
    error('No .oir files found in: %s', args.dataDir);
end

nFiles = numel(fileList);
files = {fileList.name};

% Prealloc outputs
linePeriodSeconds = nan(1, nFiles);
centroid_um       = cell(1, nFiles);
gof_all           = cell(1, nFiles);

%% --------- Read, parse, and track per file ---------
for k = 1:nFiles
    fpath = fullfile(args.dataDir, files{k});
    bf = bfopen(fpath);                        % bf{1,1}: Nx2 cell array of {plane, label}; bf{1,2}: metadata Hashtable
    planesCell = bf{1,1};
    meta       = bf{1,2};

    % Stack planes (rows=time/cycle, cols=space pixels)
    nPlanes = size(planesCell, 1);
    data = zeros(nPlanes, numel(planesCell{1,1}), 'double');
    for p = 1:nPlanes
        data(p,:) = double(planesCell{p,1});
    end

    % Try to recover line period (seconds per line). Vendors differ in keys.
    linePeriodSeconds(k) = parseLinePeriod(meta);

    % Downsample along time (rows)
    data_ds = data(1:args.downsample:end, :);
    nCycles = size(data_ds, 1);

    % Track centroid via Gaussian fit around local max on each line
    c = zeros(1, nCycles);
    r2 = nan(1, nCycles);

    W = args.roiHalfWidthPx;
    for i = 1:nCycles
        line = data_ds(i,:);
        [~, mx] = max(line);
        i1 = max(1, mx - W);
        i2 = min(numel(line), mx + W);
        seg = line(i1:i2);
        x   = (1:numel(seg))';                 %#ok<*NASGU>

        % Gaussian fit (gauss1): a1*exp(-((x-b1)/c1)^2)
        try
            [fitObj, gof] = fit((1:numel(seg))', seg', 'gauss1', 'Lower', [0, 1, 0], 'Robust', 'LAR');
            coeffs = coeffvalues(fitObj);
            mu = coeffs(2);                    % center (in segment coords)
            c(i)  = (i1 - 1) + mu;            % convert back to full line coords (pixels)
            r2(i) = gof.rsquare;
        catch
            c(i)  = NaN;
            r2(i) = NaN;
        end
    end

    % Convert to microns
    c_um = c .* args.pix2um;

    % Apply GOF filtering
    bad = r2 < args.gofThreshold;
    c_um(bad) = NaN;

    centroid_um{k} = c_um;
    gof_all{k}     = r2;
end

%% --------- Optional example preview (first file) ---------
if args.showExample
    k = 1;
    % Re-read first file planes for visualization (avoid storing full stacks for all)
    bf = bfopen(fullfile(args.dataDir, files{k}));
    planesCell = bf{1,1};
    nPlanes = size(planesCell, 1);
    A = zeros(nPlanes, numel(planesCell{1,1}), 'double');
    for p = 1:nPlanes
        A(p,:) = double(planesCell{p,1});
    end

    figure('Color','w','Name','Example line-scan preview');
    imagesc(A, [0, prctile(A(:), 99.5)]);
    colormap('gray'); axis tight; set(gca,'YDir','normal');
    ax = gca; ax.FontSize = 12;

    % Y axis: time
    if ~isnan(linePeriodSeconds(k))
        t = (0:size(A,1)-1) .* linePeriodSeconds(k);
        maxYTicks = 20;
        idx = unique(round(linspace(1, numel(t), min(maxYTicks, numel(t)))));
        yticks(idx); yticklabels(arrayfun(@(j) sprintf('%.2f', t(j)), idx, 'UniformOutput', false));
        ylabel('Time [s]', 'Interpreter','latex');
    else
        ylabel('Lines (no metadata time)', 'Interpreter','latex');
    end

    % X axis: space (microns)
    x = (0:size(A,2)-1) .* args.pix2um;
    maxXTicks = 20;
    idx = unique(round(linspace(1, numel(x), min(maxXTicks, numel(x)))));
    xticks(idx); xticklabels(arrayfun(@(j) sprintf('%.1f', x(j)), idx, 'UniformOutput', false));
    xlabel('Space [$\mu$m]', 'Interpreter','latex');
    title(sprintf('Example: %s', files{k}), 'Interpreter','none');
end

%% --------- Histogram + Gaussian fit per fiber ---------
n = nFiles;
stdDevs_um = nan(1, n);
histFits   = cell(1, n);

% Determine uniform x-limits for histograms from data range (0..p99)
allVals = [];
for k = 1:n
    v = centroid_um{k};
    v = v(~isnan(v));
    allVals = [allVals; v(:)];
end
if isempty(allVals)
    warning('All centroid values are NaN after filtering. No histograms will be shown.');
    xlims = [0, 1];
else
    hi = prctile(allVals, 99);
    xlims = [max(0, min(allVals)), max(hi, 1e-3)];
end

if n > 0
    nRows = ceil(sqrt(n));
    nCols = ceil(n / nRows);
    f = figure('Color','w','Name','Centroid histograms with Gaussian fits');
    ymax = 0;
    for k = 1:n
        v = centroid_um{k};
        v = v(~isnan(v));
        v = v(v >= xlims(1) & v <= xlims(2));
        subplot(nRows, nCols, k);
        if isempty(v)
            title(sprintf('Fiber %d (no valid data)', k));
            axis off; continue;
        end

        % Histogram
        h = histogram(v, 'BinMethod','auto'); hold on;

        % Fit Gaussian to histogram counts (bin centers vs counts)
        edges = h.BinEdges;
        xcent = edges(1:end-1) + diff(edges)/2;
        ycnt  = h.Values;
        try
            [fitObj, ~] = fit(xcent(:), ycnt(:), 'gauss1', 'Lower', [0, xlims(1), 0], 'Upper', [Inf, xlims(2), Inf], 'Robust', 'LAR');
            histFits{k} = fitObj;
            xfit = linspace(xlims(1), xlims(2), 200);
            yfit = feval(fitObj, xfit);
            plot(xfit, yfit, 'LineWidth', 1.5);
            ymax = max(ymax, max(yfit));
            % Convert Gaussian width parameter c1 to standard deviation: sigma = c1 / sqrt(2)
            stdDevs_um(k) = fitObj.c1 / sqrt(2);
        catch
            histFits{k} = [];
            stdDevs_um(k) = NaN;
        end

        title(sprintf('Fiber %d', k));
        xlabel('Centroid [\mum]');
        ylabel('Count');
        xlim(xlims);
        box on; hold off;
    end

    % harmonize y-limits
    if ymax > 0
        for k = 1:n
            subplot(nRows, nCols, k);
            yl = ylim; ylim([0, max(yl(2), ymax)]);
        end
    end
    sgtitle('Histograms with Gaussian fits');
end

%% --------- Package outputs ---------
out = struct();
out.files             = files;
out.linePeriodSeconds = linePeriodSeconds;
out.centroid_um       = centroid_um;
out.gof               = gof_all;
out.stdDevs_um        = stdDevs_um;
out.histFits          = histFits;
out.params            = args;

% Optional save
if ~isempty(args.savePrefix)
    T = table(files(:), linePeriodSeconds(:), stdDevs_um(:), ...
        'VariableNames', {'file','linePeriodSeconds','stdDev_um'});
    csvPath = fullfile(args.dataDir, sprintf('%s_summary.csv', args.savePrefix));
    matPath = fullfile(args.dataDir, sprintf('%s_results.mat', args.savePrefix));
    writetable(T, csvPath);
    save(matPath, 'out');
    fprintf('Saved: %s\nSaved: %s\n', csvPath, matPath);
end

end % main function

%% ====================== Helpers ======================

function s = withDefaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = d.(f{i});
    end
end
end

function mustExistFolder(p)
if ~isfolder(p)
    error('Folder does not exist: %s', p);
end
end

function linePeriod = parseLinePeriod(meta)
% Best-effort extraction of seconds per line from Bio-Formats metadata.
% Returns NaN if not found.
linePeriod = NaN;

try
    % Meta is typically a java.util.Hashtable
    ks = meta.keySet().toArray();
    keys = cellfun(@char, cell(ks), 'UniformOutput', false);

    candidates = [
        "speedInformation lineSpeed #2"   % legacy key observed in original code
        "speedInformation lineSpeed #1"
        "LineSpeed"
        "lineSpeed"
        "Scan Speed"
        "ScanSpeed"
        "Line Period"
        "linePeriod"
        "line period"
    ];

    valStr = "";
    for j = 1:numel(candidates)
        m = strcmpi(keys, candidates(j));
        if any(m)
            valStr = string(meta.get(candidates(j)));
            break;
        end
    end

    if strlength(valStr) == 0
        % Try any key containing "line" and "speed" or "period"
        idx = find(contains(lower(keys), "line") & (contains(lower(keys), "speed") | contains(lower(keys), "period")), 1);
        if ~isempty(idx)
            valStr = string(meta.get(keys{idx}));
        end
    end

    if strlength(valStr) > 0
        % Heuristic parsing: accept pure numbers or values with units
        % Common cases: "500.0 us/line", "0.5 ms/line", "0.0005 s", "1200 Hz" (lines per second)
        tokens = regexp(char(valStr), '([\d.]+)\s*([umk]?s|s|ms|us|Hz)?', 'tokens', 'once');
        if ~isempty(tokens)
            num = str2double(tokens{1});
            unit = '';
            if numel(tokens) >= 2 && ~isempty(tokens{2}), unit = lower(tokens{2}); end

            switch unit
                case {'s','sec','second','seconds','/line','s/line','s per line','spline'}
                    linePeriod = num;
                case {'ms','ms/line'}
                    linePeriod = num * 1e-3;
                case {'us','µs','us/line','μs','µs/line'}
                    linePeriod = num * 1e-6;
                case {'hz'} % lines per second
                    if num > 0, linePeriod = 1/num; end
                otherwise
                    % If unit missing, assume seconds if small, else microseconds if large
                    if num < 1e-1
                        linePeriod = num;
                    elseif num > 1e3
                        linePeriod = num * 1e-6;
                    else
                        linePeriod = NaN;
                    end
            end
        end
    end
catch
    linePeriod = NaN;
end

end
