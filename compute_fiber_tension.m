function Tout = compute_fiber_tension(args)
% COMPUTE_FIBER_TENSION
% Invert MSD→tension for cross-linked fibers using the odd-mode sum model:
%   msd = (2/π^4) * (L^3/ℓp) * S(φ),  with  S(φ)=Σ_{n odd} 1/(n^4 + φ n^2)
% and τ = (φ π² κ)/L², where κ = ℓp kBT.
%
% Inputs (struct 'args'):
%   args.xlsxPath     : input Excel path (required)
%   args.sheet        : sheet name or index (optional)
%   args.lp_um        : persistence length ℓp in µm       (default: 99.25347734)
%   args.kBT_pN_um    : thermal energy kBT in pN·µm       (default: 0.004114)
%   args.nMax         : largest odd mode to sum (odd up to nMax) (default: 999)
%   args.outputPath   : optional output Excel path; if empty, no file written
%   args.verbose      : true/false (default true)
%
% Expected input columns (µm and µm²):
%   - Either named exactly: "Fiber Length (um)" and "Fluctuations (um^2)"
%   - OR columns #2 and #4 will be used (for legacy spreadsheets).
%
% Output:
%   Tout : table with original columns plus:
%     "MSD_max at tau=0 (um^2)", "Feasible", "phi (dimensionless)",
%     "Tension (pN)", "Note"
%
% Usage:
%   T = compute_fiber_tension(struct( ...
%       'xlsxPath','CrossLinkedFibers_Tension_far.xlsx', ...
%       'outputPath','Fiber_Tension_Output.xlsx' ...
%   ));

%% ---- Defaults & validation ----
if nargin < 1, args = struct; end
args = withDefaults(args, struct( ...
    'xlsxPath',   '', ...
    'sheet',      [], ...
    'lp_um',      99.25347734, ...
    'kBT_pN_um',  0.004114, ...
    'nMax',       999, ...
    'outputPath', '', ...
    'verbose',    true ...
));
assert(~isempty(args.xlsxPath) && isfile(args.xlsxPath), 'xlsxPath not found.');

if args.verbose
    fprintf('Reading: %s\n', args.xlsxPath);
end
if isempty(args.sheet)
    Tin = readtable(args.xlsxPath);
else
    Tin = readtable(args.xlsxPath, 'Sheet', args.sheet);
end

%% ---- Select input columns (robust to naming) ----
varLen = "Fiber Length (um)";
varMSD = "Fluctuations (um^2)";

if all(ismember([varLen varMSD], Tin.Properties.VariableNames))
    L_um   = Tin.(varLen);
    MSD_u2 = Tin.(varMSD);
else
    % Fallback to Python script behavior: use columns #2 and #4
    assert(width(Tin) >= 4, 'Need at least 4 columns or named columns present.');
    L_um   = Tin{:, 2};
    MSD_u2 = Tin{:, 4};
    % Rename/add for clarity
    Tin.(varLen) = L_um;
    Tin.(varMSD) = MSD_u2;
end

L_um   = double(L_um(:));
MSD_u2 = double(MSD_u2(:));
N = numel(L_um);

%% ---- Constants & precompute ----
lp_um     = args.lp_um;
kBT_pN_um = args.kBT_pN_um;
odd_ns    = 1:2:args.nMax;             % 1,3,5,... up to nMax
pi2       = pi^2;
S0_limit  = (pi^4)/96;                 % sum_{odd} 1/n^4

%% ---- Per-row computation (loop for clarity & robust branching) ----
MSD_max = nan(N,1);
Feas    = false(N,1);
phi_out = nan(N,1);
tau_pN  = nan(N,1);
notes   = strings(N,1);

for i = 1:N
    L = L_um(i);
    msd = MSD_u2(i);

    if ~(L > 0) || ~(msd >= 0) || ~isfinite(L) || ~isfinite(msd)
        notes(i) = "Bad input (nonpositive or nonfinite).";
        continue;
    end

    % Zero-tension ceiling for diagnostics: msd_max = L^3/(48 ℓp)
    MSD_max(i) = (L^3) / (48*lp_um);
    Feas(i)    = (msd <= MSD_max(i) * 1.0000001);

    % Invert for φ via bracketed bisection of S(φ) - S_target = 0
    phi = phi_from_msd_scalar(L, msd, lp_um, odd_ns, S0_limit);
    phi_out(i) = phi;

    % τ = (φ π² κ) / L²,  κ = ℓp kBT (pN·µm²)
    tau = tau_from_phi_scalar(phi, L, lp_um, kBT_pN_um, pi2);
    tau_pN(i) = tau;

    if ~Feas(i)
        notes(i) = "MSD too large for any nonnegative τ (check units/data).";
    elseif ~isfinite(phi)
        notes(i) = "Numerical bracket failed (tiny MSD ⇒ huge φ).";
    else
        notes(i) = "";
    end
end

%% ---- Build output table ----
Tout = Tin;
Tout.("MSD_max at tau=0 (um^2)") = MSD_max;
Tout.("Feasible")                = Feas;
Tout.("phi (dimensionless)")     = phi_out;
Tout.("Tension (pN)")            = tau_pN;
Tout.("Note")                    = notes;

%% ---- Save if requested ----
if ~isempty(args.outputPath)
    if args.verbose, fprintf('Writing: %s\n', args.outputPath); end
    writetable(Tout, args.outputPath);
end

if args.verbose
    disp(Tout);
end

end % main

%% ================= Helpers =================
function s = withDefaults(s, d)
fn = fieldnames(d);
for k = 1:numel(fn)
    if ~isfield(s, fn{k}) || isempty(s.(fn{k}))
        s.(fn{k}) = d.(fn{k});
    end
end
end

function S = S_of_phi_scalar(phi, odd_ns)
% S(phi) = sum_{odd n} 1 / (n^4 + phi n^2)
den = odd_ns.^2;
S = sum(1 ./ (den.^2 + phi * den));
end

function phi = phi_from_msd_scalar(L_um, msd_um2, lp_um, odd_ns, S0_limit)
% Solve S(phi) = π^4 ℓp msd / (2 L^3) for φ >= 0 via bracketed bisection.
if ~(L_um > 0) || ~(msd_um2 >= 0) || ~isfinite(L_um) || ~isfinite(msd_um2)
    phi = NaN; return;
end

S_target = (pi^4) * lp_um * msd_um2 / (2 * (L_um^3));

% If target exceeds S(0), no physical nonnegative tension exists.
if S_target > S0_limit * 1.0000001
    phi = NaN; return;
end

f = @(phi) S_of_phi_scalar(phi, odd_ns) - S_target;

lo = 0.0; flo = f(lo);
if ~isfinite(flo), phi = NaN; return; end

hi = 1.0; fhi = f(hi);
tries = 0;
while fhi > 0.0 && tries < 60
    hi = hi * 2.0;
    fhi = f(hi);
    tries = tries + 1;
end

if fhi > 0.0
    % Could not bracket (msd→0 ⇒ φ→∞). Report Inf (caller will note).
    phi = Inf; return;
end

% Bisection
for iter = 1:200
    mid = 0.5*(lo + hi);
    fm  = f(mid);
    if abs(fm) < 1e-13
        phi = mid; return;
    end
    if fm > 0.0
        lo = mid;
    else
        hi = mid;
    end
    if (hi - lo) < 1e-12 * (1.0 + lo + hi)
        phi = 0.5*(lo + hi); return;
    end
end

phi = 0.5*(lo + hi);
end

function tau = tau_from_phi_scalar(phi, L_um, lp_um, kBT_pN_um, pi2)
% τ [pN] = (φ π² κ)/L², κ = ℓp kBT [pN·µm²]
if ~isfinite(phi) || ~(L_um > 0)
    tau = NaN; return;
end
kappa = lp_um * kBT_pN_um;
tau = (phi * pi2 * kappa) / (L_um^2);
end
