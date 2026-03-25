%% AOC Multiverse — Subject-Level Data Preparation (Sternberg & N-Back)
% Computes subject-level (trial-averaged) alpha and gaze features for all
% multiverse dimensions. Directly comparable to standard subject-level
% analyses (rainclouds, etc.).
%
% Key differences from trial-level (AOC_multiverse_prep.m):
%   - EEG power spectra averaged across trials BEFORE alpha extraction
%   - FOOOF fitted to trial-averaged spectra (not per-trial → ~75x faster)
%   - EEG baseline correction applied to averaged spectra
%   - Gaze metrics computed per-trial, then averaged to subject means
%   - Output: 1 row per subject × condition × universe (~360k rows vs ~31M)
%
% Writes multiverse_sternberg_subject.csv and multiverse_nback_subject.csv.
%
% Decision grid (7 dimensions, 1440 universes per task):
%   Electrodes:     posterior, occipital (2)
%   1/f:            FOOOFed, non-FOOOFed (2)
%   Latency:        0-500, 0-1000, 0-2000, 1000-2000 ms (4)
%   Alpha band:     canonical 8-14 Hz, IAF (2)
%   Gaze measure:   SPL, velocity, microsaccades, BCEA, gaze deviation (5)
%   EEG baseline:   raw, dB, pct_change [-0.5 -0.25] s (3)
%   Gaze baseline:  raw, dB, pct_change [-0.5 -0.25] s (3)
%
% Run AFTER preprocessing (requires the same feature files as trial-level).

disp(upper('=== AOC MULTIVERSE SUBJECT-LEVEL PREP START ==='))

%% Setup: run startup and setup FIRST (startup may clear the workspace)
path_preproc = [];
subjects = {};
try
    if exist('startup', 'file')
        startup
        disp(upper('Startup run.'))
    end
    if exist('setup', 'file')
        [subjects, paths, ~, ~] = setup('AOC');
        path_preproc = paths.features;
        disp(upper(['Setup: ' num2str(length(subjects)) ' subjects, path_preproc = ' path_preproc]))
    end
catch ME
    disp(upper(['Setup failed: ' ME.message ' — using base_features for subject list.']))
end

%% Paths (Science Cloud vs Mac) — defined AFTER startup/setup so clear all cannot wipe them
if ispc
    base_data_in = 'W:\Students\Arne\AOC';
    base_data_out = 'W:\Students\Arne\MAG';
    base_features = fullfile(base_data_in, 'data', 'features');
else
    base_data_in = '/Volumes/g_psyplafor_methlab$/Students/Arne/AOC';
    base_data_out = '/Volumes/g_psyplafor_methlab$/Students/Arne/MAG';
    base_features = fullfile(base_data_in, 'data', 'features');
end
out_dir = fullfile(base_data_out, 'data', 'multiverse');
if ~isfolder(out_dir), mkdir(out_dir); end
disp(upper(['Paths: base_data_in = ' base_data_in ', base_data_out = ' base_data_out ', out_dir (CSV) = ' out_dir]))

if isempty(subjects)
    dirs = dir(base_features);
    dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {dirs.name};
    disp(upper(['Subject list from base_features: ' num2str(length(subjects)) ' folders.']))
end
if isempty(path_preproc)
    path_preproc = base_features;
    disp(upper('Using base_features as path_preproc.'))
end

%% Decision options (7 dimensions)
electrodes_opts    = {'posterior', 'occipital'};
fooof_opts         = {'FOOOFed', 'nonFOOOFed'};
latency_opts       = {'0_500ms', '0_1000ms', '0_2000ms', '1000_2000ms'};
alpha_opts         = {'canonical', 'IAF'};
gaze_opts          = {'scan_path_length', 'gaze_velocity', 'microsaccades', 'BCEA', 'gaze_deviation'};
gaze_col_map       = [2, 3, 4, 5, 1];  % column indices into [gaze_dev, spl, vel, ms, bcea]
baseline_eeg_opts  = {'raw', 'dB', 'pct_change'};
baseline_gaze_opts = {'raw', 'dB', 'pct_change'};
n_elec = 2; n_fooof = 2; n_lat = 4; n_alpha = 2; n_gaze = 5;
n_bl_eeg = 3; n_bl_gaze = 3;
n_universes = n_elec * n_fooof * n_lat * n_alpha * n_gaze * n_bl_eeg * n_bl_gaze;
alphaRange = [8 14];
disp(upper(['Decision grid: ' num2str(n_universes) ' universes per task (7 dimensions).']))

%% Get channel labels and indices
disp(upper('Resolving channel sets (posterior, occipital)...'))
labels_master = [];
for s = 1:length(subjects)
    sid = subjects{s};
    eeg_dir = fullfile(path_preproc, sid, 'eeg');
    early_file = fullfile(eeg_dir, 'power_stern_early_trials.mat');
    if ~isfile(early_file)
        early_file = fullfile(base_features, sid, 'eeg', 'power_stern_early_trials.mat');
    end
    if isfile(early_file)
        load(early_file, 'powload2_early');
        labels_master = powload2_early.label;
        break
    end
end
if isempty(labels_master)
    error('Could not load power_stern_early_trials from any subject.')
end
[idx_post, idx_occ] = get_channel_indices(labels_master);
ch_sets = {idx_post, idx_occ};
ch_label_sets = {labels_master(idx_post), labels_master(idx_occ)};
disp(upper(['Channels: post=' num2str(length(idx_post)) ...
  ' occ=' num2str(length(idx_occ))]))

%% FOOOF R² output directory
r2_dir = fullfile(base_data_out, 'data', 'controls', 'multiverse');
if ~isfolder(r2_dir), mkdir(r2_dir); end

%% Build Sternberg multiverse table
disp(upper('--- STERNBERG TASK ---'))
[tbl_s, r2_s] = build_task_multiverse_subject('sternberg', subjects, path_preproc, base_features, ...
    ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
    baseline_eeg_opts, baseline_gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange);
disp(upper(['Sternberg table rows: ' num2str(height(tbl_s))]))
if height(tbl_s) == 0
  error('AOC_multiverse_prep_subject:NoData', 'Sternberg table is empty.')
end
writetable(tbl_s, fullfile(out_dir, 'multiverse_sternberg_subject.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_sternberg_subject.csv')]))
if height(r2_s) > 0
  writetable(r2_s, fullfile(r2_dir, 'fooof_r2_sternberg_subject.csv'));
  disp(upper(['Written: ' fullfile(r2_dir, 'fooof_r2_sternberg_subject.csv') ' (' num2str(height(r2_s)) ' FOOOF fits)']))
end

%% Build N-back multiverse table
disp(upper('--- N-BACK TASK ---'))
[tbl_n, r2_n] = build_task_multiverse_subject('nback', subjects, path_preproc, base_features, ...
    ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
    baseline_eeg_opts, baseline_gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange);
disp(upper(['N-back table rows: ' num2str(height(tbl_n))]))
if height(tbl_n) == 0
  error('AOC_multiverse_prep_subject:NoData', 'N-back table is empty.')
end
writetable(tbl_n, fullfile(out_dir, 'multiverse_nback_subject.csv'));
disp(upper(['Written: ' fullfile(out_dir, 'multiverse_nback_subject.csv')]))
if height(r2_n) > 0
  writetable(r2_n, fullfile(r2_dir, 'fooof_r2_nback_subject.csv'));
  disp(upper(['Written: ' fullfile(r2_dir, 'fooof_r2_nback_subject.csv') ' (' num2str(height(r2_n)) ' FOOOF fits)']))
end

disp(upper('=== AOC MULTIVERSE SUBJECT-LEVEL PREP DONE ==='))

%% ========== LOCAL FUNCTIONS ==========

function [idx_post, idx_occ] = get_channel_indices(labels)
  n = length(labels);
  idx_occ = []; idx_pari = [];
  for i = 1:n
    lb = labels{i};
    if contains(lb, 'O') || contains(lb, 'I'), idx_occ(end+1) = i; end
    if contains(lb, 'P') && ~contains(lb, 'F'), idx_pari(end+1) = i; end
  end
  idx_post = unique([idx_occ, idx_pari]);
  if isempty(idx_post), idx_post = idx_occ; end
  disp(['    Posterior (' num2str(length(idx_post)) '): ' strjoin(labels(idx_post), ', ')])
  disp(['    Occipital (' num2str(length(idx_occ)) '): ' strjoin(labels(idx_occ), ', ')])
end

function a = bandpower_from_spectrum(mean_spectrum, freq, chIdx, band)
  % Extract mean alpha power from a trial-averaged spectrum (chan × freq)
  if isempty(chIdx), a = NaN; return, end
  f = freq(:);
  bandIdx = f >= band(1) & f <= band(2);
  if sum(bandIdx) == 0, a = NaN; return, end
  x = mean_spectrum(chIdx, bandIdx);
  a = mean(x(:), 'omitnan');
end

function IAF_band = get_IAF_band(pow_full, chIdx, alphaRange)
  if isempty(chIdx), IAF_band = alphaRange; return, end
  spec = squeeze(mean(mean(pow_full.powspctrm(:, chIdx, :), 1), 2));
  f = pow_full.freq(:);
  aIdx = find(f >= alphaRange(1) & f <= alphaRange(2));
  alphaF = f(aIdx); alphaP = spec(aIdx);
  [pks, locs] = findpeaks(double(alphaP));
  if isempty(pks), IAF_band = alphaRange; return, end
  [~, im] = max(pks);
  IAF = alphaF(locs(im));
  if IAF <= alphaRange(1) || IAF >= alphaRange(2), IAF_band = alphaRange; return, end
  IAF_band = [max(f(1), IAF-4), min(f(end), IAF+2)];
end

function bcea = compute_bcea(xw, yw)
  k95 = chi2inv(0.95, 2) / 2;
  xw = double(xw(:)); yw = double(yw(:));
  if numel(xw) < 10, bcea = NaN; return, end
  sx = std(xw); sy = std(yw);
  if sx == 0 || sy == 0, bcea = NaN; return, end
  rho = corr(xw, yw);
  bcea = 2 * k95 * pi * sx * sy * sqrt(1 - rho^2);
end

function [gaze_raw, gaze_db, gaze_pct] = compute_gaze_one_window(x, y, t, tw, dur, fsample, t_base)
  gaze_raw = nan(1, 5); gaze_db = nan(1, 5); gaze_pct = nan(1, 5);
  idx = t >= tw(1) & t <= tw(2);
  xw = x(idx); yw = y(idx);
  validxy = isfinite(xw) & isfinite(yw);
  xw = xw(validxy); yw = yw(validxy);
  if length(xw) < 50, return, end
  gaze_dev = mean(sqrt((xw - 400).^2 + (yw - 300).^2), 'omitnan');
  dx = diff(xw); dy = diff(yw);
  spl = sum(sqrt(dx.^2 + dy.^2), 'omitnan');
  vel = spl / dur;
  try
    [ms_rate, ms_det] = detect_microsaccades(fsample, [xw; yw], length(xw));
    if exist('ms_det', 'var') && isfield(ms_det, 'Onset')
      ms_cnt = numel(ms_det.Onset);
    else
      ms_cnt = ms_rate * dur;
    end
  catch, ms_cnt = NaN;
  end
  bcea_val = compute_bcea(xw, yw);
  gaze_raw = [gaze_dev, spl, vel, ms_cnt, bcea_val];
  idx_b = t >= t_base(1) & t <= t_base(2);
  xb = x(idx_b); yb = y(idx_b);
  validb = isfinite(xb) & isfinite(yb);
  xb = xb(validb); yb = yb(validb);
  dur_base = t_base(2) - t_base(1);
  if length(xb) < 20, return, end
  gaze_dev_b = mean(sqrt((xb - 400).^2 + (yb - 300).^2), 'omitnan');
  dx_b = diff(xb); dy_b = diff(yb);
  spl_b = sum(sqrt(dx_b.^2 + dy_b.^2), 'omitnan');
  vel_b = spl_b / dur_base;
  try
    [ms_rate_b, ms_det_b] = detect_microsaccades(fsample, [xb; yb], length(xb));
    if exist('ms_det_b', 'var') && isfield(ms_det_b, 'Onset')
      ms_cnt_b = numel(ms_det_b.Onset);
    else
      ms_cnt_b = ms_rate_b * dur_base;
    end
  catch, ms_cnt_b = NaN;
  end
  bcea_b = compute_bcea(xb, yb);
  base_vals = [gaze_dev_b, spl_b, vel_b, ms_cnt_b, bcea_b];
  for g = 1:5
    if isfinite(base_vals(g)) && base_vals(g) > 0 && isfinite(gaze_raw(g))
      gaze_pct(g) = (gaze_raw(g) - base_vals(g)) / abs(base_vals(g)) * 100;
      if gaze_raw(g) > 0
        gaze_db(g) = 10 * log10(gaze_raw(g) / base_vals(g));
      end
    end
  end
end

function pow = get_power_from_TFR(tfr_all, ind_cond, latency)
  cfg = []; cfg.latency = latency; cfg.trials = ind_cond;
  sel = ft_selectdata(cfg, tfr_all);
  sel.powspctrm = mean(sel.powspctrm, 4);
  if isfield(sel, 'time'), sel = rmfield(sel, 'time'); end
  sel.dimord = 'rpt_chan_freq';
  pow = sel;
end

function pow = compute_pow_cond_window(data_td, cond_inds, time_win, cfg_freq)
  pow = [];
  if isempty(cond_inds), return, end
  try
    cfg_sel = []; cfg_sel.trials = cond_inds;
    d = ft_selectdata(cfg_sel, data_td);
    cfg_lat = []; cfg_lat.latency = time_win;
    d = ft_selectdata(cfg_lat, d);
    pow = ft_freqanalysis(cfg_freq, d);
  catch ME
    disp(upper(['    WARNING: compute_pow_cond_window failed: ' ME.message]))
    pow = [];
  end
end

function [f_osc, f_axis, ap_offset, ap_exponent, fooof_r2, fooof_err] = run_fooof_subject(data_td, cond_inds, ch_labels, time_win, cfg_fooof)
  % Run FOOOF on trial-averaged spectrum (all condition trials at once).
  % cfg_fooof.keeptrials = 'no' → ft_freqanalysis internally averages trials.
  % Also returns aperiodic parameters (offset, exponent) and fit quality (R², error).
  f_osc = []; f_axis = []; ap_offset = NaN; ap_exponent = NaN; fooof_r2 = NaN; fooof_err = NaN;
  if isempty(data_td) || isempty(ch_labels) || isempty(cond_inds), return, end
  try
    cfg_sel = []; cfg_sel.trials = cond_inds;
    d = ft_selectdata(cfg_sel, data_td);
    cfg_lat = []; cfg_lat.latency = time_win;
    d = ft_selectdata(cfg_lat, d);
    cfg_ch = []; cfg_ch.channel = ch_labels; cfg_ch.avgoverchan = 'yes';
    d = ft_selectdata(cfg_ch, d);
    d.label = {'ROI'};
    if ~exist('ft_freqanalysis_Arne_FOOOF', 'file'), return, end
    out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, d);
    f_axis = out.freq(:)';
    if iscell(out.fooofparams), rep = out.fooofparams{1}; else, rep = out.fooofparams; end
    if isfield(rep, 'aperiodic_params') && numel(rep.aperiodic_params) >= 2
      ap_offset = rep.aperiodic_params(1);
      ap_exponent = rep.aperiodic_params(2);
    end
    if isfield(rep, 'r_squared'), fooof_r2 = rep.r_squared; end
    if isfield(rep, 'error'), fooof_err = rep.error; end
    if isfield(rep, 'fooofed_spectrum') && ~isempty(rep.fooofed_spectrum)
      f_osc = rep.fooofed_spectrum(:)';
    elseif isfield(rep, 'peak_params') && ~isempty(rep.peak_params)
      pk = rep.peak_params;
      g = zeros(size(out.freq(:)'));
      for p = 1:size(pk,1)
        g = g + pk(p,2) .* exp(-(out.freq(:)' - pk(p,1)).^2 ./ (2*pk(p,3)^2));
      end
      f_osc = g;
    end
  catch ME
    disp(upper(['    WARNING: run_fooof_subject failed: ' ME.message]))
    f_osc = []; f_axis = [];
  end
end

function mean_spec = average_spectrum(pow)
  % Average trial-level spectra → 1 × chan × freq (handles NaN)
  mean_spec = [];
  if isempty(pow) || ~isfield(pow, 'powspctrm'), return, end
  s = pow;
  s.powspctrm = mean(pow.powspctrm, 1, 'omitnan');  % 1 × chan × freq
  s.dimord = 'rpt_chan_freq';
  mean_spec = s;
end

function pow_bl = compute_db_baseline(mean_task, mean_base)
  % dB baseline on averaged spectra: 10*log10(task_mean / base_mean)
  pow_bl = [];
  if isempty(mean_task) || isempty(mean_base), return, end
  pow_bl = mean_task;
  base = mean_base.powspctrm;  % 1 × chan × freq
  if length(mean_base.freq) ~= length(mean_task.freq) || ...
      max(abs(mean_base.freq - mean_task.freq)) > 1e-6
    disp(upper('    NOTE: Freq axis mismatch. Interpolating baseline.'))
    base_aligned = nan(1, size(base, 2), length(mean_task.freq));
    for ch = 1:size(base, 2)
      base_aligned(1, ch, :) = interp1(mean_base.freq, squeeze(base(1, ch, :)), ...
          mean_task.freq, 'linear', NaN);
    end
    base = base_aligned;
  end
  base(base <= 0) = NaN;
  pow_bl.powspctrm = 10 * log10(bsxfun(@rdivide, mean_task.powspctrm, base));
end

function pow_bl = compute_pct_baseline(mean_task, mean_base)
  % Percentage-change baseline on averaged spectra
  pow_bl = [];
  if isempty(mean_task) || isempty(mean_base), return, end
  pow_bl = mean_task;
  base = mean_base.powspctrm;
  if length(mean_base.freq) ~= length(mean_task.freq) || ...
      max(abs(mean_base.freq - mean_task.freq)) > 1e-6
    disp(upper('    NOTE: Freq axis mismatch. Interpolating baseline.'))
    base_aligned = nan(1, size(base, 2), length(mean_task.freq));
    for ch = 1:size(base, 2)
      base_aligned(1, ch, :) = interp1(mean_base.freq, squeeze(base(1, ch, :)), ...
          mean_task.freq, 'linear', NaN);
    end
    base = base_aligned;
  end
  base(base <= 0) = NaN;
  pow_bl.powspctrm = bsxfun(@rdivide, ...
      bsxfun(@minus, mean_task.powspctrm, base), abs(base)) * 100;
end

%% ========== MAIN BUILD FUNCTION ==========

function [tbl, r2_tbl] = build_task_multiverse_subject(task_name, subjects, path_preproc, base_features, ...
    ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
    baseline_eeg_opts, baseline_gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange)

  disp(upper(['Building SUBJECT-LEVEL multiverse table for task: ' task_name ' (' num2str(n_universes) ' universes)']))
  if strcmp(task_name, 'sternberg')
    cond_codes = [22 24 26]; cond_vals = [2 4 6];
    early_file = 'power_stern_early_trials.mat'; full_file = 'power_stern_full_trials.mat';
    tfr_file = 'tfr_stern_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_sternberg.mat';
    eeg_td_file = 'dataEEG_sternberg.mat';
    et_file = 'dataET_sternberg.mat'; et_var = 'dataETlong';
  else
    cond_codes = [21 22 23]; cond_vals = [1 2 3];
    early_file = 'power_nback_early_trials.mat'; full_file = 'power_nback_full_trials.mat';
    tfr_file = 'tfr_nback_trials.mat'; eeg_tfr_file = 'dataEEG_TFR_nback.mat';
    eeg_td_file = 'dataEEG_nback.mat';
    et_file = 'dataET_nback.mat'; et_var = 'dataETlong';
  end

  cfg_hann = []; cfg_hann.method = 'mtmfft'; cfg_hann.taper = 'hanning';
  cfg_hann.foilim = [2 40]; cfg_hann.pad = 5; cfg_hann.keeptrials = 'yes';

  cfg_fooof = []; cfg_fooof.method = 'mtmfft'; cfg_fooof.taper = 'hanning';
  cfg_fooof.foilim = [2 40]; cfg_fooof.pad = 5; cfg_fooof.output = 'fooof'; cfg_fooof.keeptrials = 'no';

  lat_windows = {[0 0.5], [0 1], [0 2], [1 2]};
  lat_durs    = [0.5, 1, 2, 1];
  t_base = [-0.5 -0.25];
  fsample = 500;

  % Preallocate output arrays (no Trial column — 1 row per subject × condition × universe)
  est_rows = length(subjects) * 3 * n_universes;
  task_cell = cell(est_rows, 1); u_id_cell = nan(est_rows, 1);
  elec_cell = cell(est_rows, 1); fooof_cell = cell(est_rows, 1); lat_cell = cell(est_rows, 1);
  alpha_type_cell = cell(est_rows, 1); gaze_meas_cell = cell(est_rows, 1);
  bl_eeg_cell = cell(est_rows, 1); bl_gaze_cell = cell(est_rows, 1);
  subjectID_cell = nan(est_rows, 1); Condition_cell = nan(est_rows, 1);
  alpha_val_cell = nan(est_rows, 1); gaze_val_cell = nan(est_rows, 1);
  n_trials_cell = nan(est_rows, 1);
  ap_off_val_cell = nan(est_rows, 1); ap_exp_val_cell = nan(est_rows, 1);
  row_idx = 0;

  % FOOOF R² collection (compact: one row per FOOOF call)
  r2_sid = []; r2_cond = []; r2_lat = {}; r2_elec = {};
  r2_val = []; r2_err = []; r2_exp = []; r2_off = [];

  for s = 1:length(subjects)
    sid = subjects{s};
    sid_num = str2double(sid);
    if isnan(sid_num), continue; end
    disp(upper(['Subject ' sid ' (' num2str(s) '/' num2str(length(subjects)) ')']))

    eeg_dir = fullfile(path_preproc, sid, 'eeg');
    gaze_dir = fullfile(base_features, sid, 'gaze');
    if ~isfolder(eeg_dir), eeg_dir = fullfile(base_features, sid, 'eeg'); end
    if ~isfolder(gaze_dir), gaze_dir = fullfile(path_preproc, sid, 'gaze'); end

    %% ====== EEG: Load pre-computed trial-level power (Hanning, raw) ======
    disp(upper(['  Loading EEG: ' eeg_dir]))
    early_path = fullfile(eeg_dir, early_file);
    full_path  = fullfile(eeg_dir, full_file);
    if ~isfile(early_path) || ~isfile(full_path)
      disp(upper('  Skip: missing power files.')); continue
    end
    load(early_path); load(full_path);
    disp(upper('  Loaded power early (0-1 s) and full (0-2 s).'))

    if strcmp(task_name, 'sternberg')
      pow_early = {powload2_early, powload4_early, powload6_early};
      pow_full  = {powload2_full,  powload4_full,  powload6_full};
    else
      pow_early = {powload1_early, powload2_early, powload3_early};
      pow_full  = {powload1_full,  powload2_full,  powload3_full};
    end

    % Validate pow_full: N-back TFR uses cfg.toi up to 2.25s, but the 1s
    % convolution window causes edge NaN at t>=1.75s. mean() over time then
    % propagates NaN to the entire trial-averaged spectrum. If all NaN,
    % clear so the fallback (direct FFT on time-domain data) kicks in.
    for c = 1:length(pow_full)
      if ~isempty(pow_full{c}) && isfield(pow_full{c}, 'powspctrm') && ...
          all(isnan(pow_full{c}.powspctrm(:)))
        disp(upper(['  WARNING: pow_full{' num2str(c) '} is all NaN (TFR edge effect). Will recompute from time-domain.']))
        pow_full{c} = [];
      end
    end

    %% ====== Load TFR ======
    tfr_path = fullfile(eeg_dir, tfr_file);
    eeg_tfr_path = fullfile(eeg_dir, eeg_tfr_file);
    eeg_td_path = fullfile(eeg_dir, eeg_td_file);

    tfr_loaded = [];
    if isfile(tfr_path)
      disp(upper('  Loading precomputed TFR (Hanning).'))
      load(tfr_path, 'tfr_all');
      tfr_loaded = tfr_all;
    end

    %% ====== Load time-domain EEG (for Hanning fallback + FOOOF) ======
    data_td = [];
    if isfile(eeg_tfr_path)
      disp(upper('  Loading time-domain EEG (dataEEG_TFR).'))
      load(eeg_tfr_path, 'dataTFR');
      data_td = dataTFR;
    elseif isfile(eeg_td_path)
      disp(upper('  Loading time-domain EEG (dataEEG).'))
      load(eeg_td_path, 'dataEEG');
      if ~exist('dataEEG', 'var')
        load(eeg_td_path, 'dataTFR'); data_td = dataTFR;
      else
        data_td = dataEEG;
      end
    end
    if isempty(data_td) && isempty(tfr_loaded)
      disp(upper('  WARNING: No TFR or time-domain EEG. Some latencies and all FOOOF will be NaN.'))
    end

    %% ====== Subject IAF (from Hanning full or early, occipital, raw) ======
    iaf_pow = pow_full{1};
    if isempty(iaf_pow), iaf_pow = pow_early{1}; end
    if isempty(iaf_pow)
      disp(upper('  WARNING: No valid power data for IAF. Using default alpha range.'))
      IAF_band = alphaRange;
    else
      disp(upper('  Computing subject IAF from occipital power.'))
      IAF_band = get_IAF_band(iaf_pow, ch_sets{2}, alphaRange);
    end

    %% ====== Load gaze ======
    disp(upper(['  Loading gaze: ' gaze_dir]))
    et_path = fullfile(gaze_dir, et_file);
    if ~isfile(et_path)
      disp(upper('  Skip: missing ET file.')); continue
    end
    load(et_path);
    if ~exist(et_var, 'var'), et_var = 'dataET'; end
    if ~exist(et_var, 'var'), disp(upper('  Skip: ET variable not found.')); continue; end
    dataET = eval(et_var);
    disp(upper(['  Gaze loaded: ' num2str(size(dataET.trialinfo, 1)) ' trials.']))

    trial_num = dataET.trialinfo(:,2);
    cond_code = dataET.trialinfo(:,1);
    if strcmp(task_name, 'sternberg'), cond_code = cond_code - 20; end
    if strcmp(task_name, 'nback'),     cond_code = cond_code - 20; end

    %% ====== Loop over conditions ======
    for cond = 1:3
      cval = cond_vals(cond);
      trl_idx = find(cond_code == cval);
      if isempty(trl_idx), continue; end

      n_trials_cond = length(trl_idx);
      n_eeg_trials = size(pow_early{cond}.powspctrm, 1);
      if n_trials_cond ~= n_eeg_trials
        disp(upper(['  WARNING: ET=' num2str(n_trials_cond) ' vs EEG=' num2str(n_eeg_trials) ...
            ' for cond ' num2str(cval) '. Using min.']))
        n_trials_cond = min(n_trials_cond, n_eeg_trials);
        trl_idx = trl_idx(1:n_trials_cond);
      end
      disp(upper(['  Condition ' num2str(cval) ': ' num2str(n_trials_cond) ' trials.']))

      %% ====== 1. Build trial-level power structs, then AVERAGE across trials ======
      pow_s = cell(4, 1);  % trial-level raw power per latency
      pow_s{2} = pow_early{cond};  % 0-1s
      pow_s{3} = pow_full{cond};   % 0-2s

      if ~isempty(tfr_loaded)
        cond_inds_tfr = find(tfr_loaded.trialinfo(:,1) == cond_codes(cond));
        if isempty(pow_s{1})
          pow_s{1} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{1});
        end
        if isempty(pow_s{4})
          pow_s{4} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{4});
        end
      end

      cond_inds_td = [];
      if ~isempty(data_td)
        cond_inds_td = find(data_td.trialinfo(:,1) == cond_codes(cond));
        for il = 1:4
          if isempty(pow_s{il})
            disp(upper(['    Computing power from time-domain: lat=' num2str(il)]))
            pow_s{il} = compute_pow_cond_window(data_td, cond_inds_td, lat_windows{il}, cfg_hann);
          end
        end
      end

      % Average power across trials → subject-mean spectra
      disp(upper('    Averaging power spectra across trials.'))
      mean_raw = cell(4, 1);
      for il = 1:4
        mean_raw{il} = average_spectrum(pow_s{il});
      end

      % Baseline spectrum (averaged)
      pow_base_trials = [];
      if ~isempty(tfr_loaded)
        try
          pow_base_trials = get_power_from_TFR(tfr_loaded, cond_inds_tfr, t_base);
        catch
          pow_base_trials = [];
        end
      end
      if isempty(pow_base_trials) && ~isempty(data_td) && ~isempty(cond_inds_td)
        pow_base_trials = compute_pow_cond_window(data_td, cond_inds_td, t_base, cfg_hann);
      end
      mean_base = average_spectrum(pow_base_trials);

      % Baseline-correct averaged spectra
      mean_db  = cell(4, 1);
      mean_pct = cell(4, 1);
      for il = 1:4
        if ~isempty(mean_raw{il}) && ~isempty(mean_base)
          mean_db{il}  = compute_db_baseline(mean_raw{il}, mean_base);
          mean_pct{il} = compute_pct_baseline(mean_raw{il}, mean_base);
        end
      end

      %% ====== 2. Extract non-FOOOF alpha from averaged spectra ======
      % alpha_nf{lat, bl}(elec_alpha_col) — scalar values
      n_cols = n_elec * n_alpha;
      alpha_nf = cell(n_lat, n_bl_eeg);
      for il = 1:n_lat; for ibl = 1:n_bl_eeg
        alpha_nf{il, ibl} = nan(1, n_cols);
      end; end

      for il = 1:n_lat
        for ibl = 1:n_bl_eeg
          if ibl == 1, p = mean_raw{il}; elseif ibl == 2, p = mean_db{il}; else, p = mean_pct{il}; end
          if isempty(p), continue; end
          spec = squeeze(p.powspctrm(1, :, :));  % chan × freq
          for ie = 1:n_elec
            for ia = 1:n_alpha
              col = (ie-1)*n_alpha + ia;
              if ia == 1, band = alphaRange; else, band = IAF_band; end
              alpha_nf{il, ibl}(col) = bandpower_from_spectrum(spec, p.freq, ch_sets{ie}, band);
            end
          end
        end
      end

      %% ====== 3. FOOOF on trial-averaged data (1 call per latency × electrode) ======
      % fooof_alpha{lat}(elec_alpha_col) — scalar values
      fooof_alpha = cell(n_lat, 1);
      ap_offset_subj = cell(n_lat, 1);
      ap_exp_subj = cell(n_lat, 1);
      for il = 1:n_lat
        fooof_alpha{il} = nan(1, n_cols);
        ap_offset_subj{il} = nan(1, n_elec);
        ap_exp_subj{il} = nan(1, n_elec);
      end

      if ~isempty(data_td) && ~isempty(cond_inds_td)
        disp(upper('    Running FOOOF on trial-averaged spectra.'))
        for il = 1:n_lat
          for ie = 1:n_elec
            [f_osc, f_axis, ap_off, ap_exp, f_r2, f_err] = run_fooof_subject(data_td, cond_inds_td, ...
                ch_label_sets{ie}, lat_windows{il}, cfg_fooof);
            if isempty(f_osc)
              ap_offset_subj{il}(ie) = ap_off;
              ap_exp_subj{il}(ie) = ap_exp;
              continue
            end
            ap_offset_subj{il}(ie) = ap_off;
            ap_exp_subj{il}(ie) = ap_exp;
            r2_sid(end+1,1) = sid_num; r2_cond(end+1,1) = cval;
            r2_lat{end+1,1} = latency_opts{il}; r2_elec{end+1,1} = electrodes_opts{ie};
            r2_val(end+1,1) = f_r2; r2_err(end+1,1) = f_err;
            r2_exp(end+1,1) = ap_exp; r2_off(end+1,1) = ap_off;
            for ia = 1:n_alpha
              col = (ie-1)*n_alpha + ia;
              if ia == 1, band = alphaRange; else, band = IAF_band; end
              bandIdx = f_axis >= band(1) & f_axis <= band(2);
              if sum(bandIdx) > 0
                fooof_alpha{il}(col) = mean(f_osc(bandIdx), 'omitnan');
              end
            end
          end
        end
      else
        disp(upper('    WARNING: No time-domain EEG → all FOOOFed alpha will be NaN.'))
      end

      %% ====== 4. Gaze: per-trial computation, then average to subject mean ======
      disp(upper('    Computing per-trial gaze metrics and averaging.'))
      gaze_raw_mean = cell(n_lat, 1);
      gaze_db_mean  = cell(n_lat, 1);
      gaze_pct_mean = cell(n_lat, 1);

      for il = 1:n_lat
        gaze_raw_tr = nan(n_trials_cond, 5);
        gaze_db_tr  = nan(n_trials_cond, 5);
        gaze_pct_tr = nan(n_trials_cond, 5);
        for tr = 1:n_trials_cond
          trl_glob = trl_idx(tr);
          raw_et = dataET.trial{trl_glob};
          t = dataET.time{trl_glob};
          if size(raw_et,1) < 2, continue; end
          x = raw_et(1,:); y = raw_et(2,:);
          if length(t) ~= length(x), continue; end
          [gaze_raw_tr(tr,:), gaze_db_tr(tr,:), gaze_pct_tr(tr,:)] = ...
            compute_gaze_one_window(x, y, t, lat_windows{il}, lat_durs(il), fsample, t_base);
        end
        gaze_raw_mean{il} = mean(gaze_raw_tr, 1, 'omitnan');
        gaze_db_mean{il}  = mean(gaze_db_tr,  1, 'omitnan');
        gaze_pct_mean{il} = mean(gaze_pct_tr, 1, 'omitnan');
      end

      %% ====== 5. Append 1 row per universe ======
      disp(upper(['    Appending ' num2str(n_universes) ' rows for subject ' sid ' cond ' num2str(cval) '.']))
      for u = 1:n_universes
        [ie, ifo, il, ia, ig, ibe, ibg] = ind2sub( ...
          [n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze], u);

        row_idx = row_idx + 1;

        task_cell{row_idx}       = task_name;
        u_id_cell(row_idx)       = u;
        elec_cell{row_idx}       = electrodes_opts{ie};
        fooof_cell{row_idx}      = fooof_opts{ifo};
        lat_cell{row_idx}        = latency_opts{il};
        alpha_type_cell{row_idx} = alpha_opts{ia};
        gaze_meas_cell{row_idx}  = gaze_opts{ig};
        bl_eeg_cell{row_idx}     = baseline_eeg_opts{ibe};
        bl_gaze_cell{row_idx}    = baseline_gaze_opts{ibg};
        subjectID_cell(row_idx)  = sid_num;
        Condition_cell(row_idx)  = cval;
        n_trials_cell(row_idx)   = n_trials_cond;

        % Alpha: FOOOF is baseline-independent → same for all EEG baselines
        col_alpha = (ie-1)*n_alpha + ia;
        if ifo == 1  % FOOOFed
          alpha_val_cell(row_idx) = fooof_alpha{il}(col_alpha);
        else         % nonFOOOFed
          alpha_val_cell(row_idx) = alpha_nf{il, ibe}(col_alpha);
        end

        % Gaze lookup
        gcol = gaze_col_map(ig);
        if ibg == 1
          gaze_val_cell(row_idx) = gaze_raw_mean{il}(gcol);
        elseif ibg == 2
          gaze_val_cell(row_idx) = gaze_db_mean{il}(gcol);
        else
          gaze_val_cell(row_idx) = gaze_pct_mean{il}(gcol);
        end

        % Aperiodic lookup (only valid for FOOOFed universes)
        if ifo == 1
          ap_off_val_cell(row_idx) = ap_offset_subj{il}(ie);
          ap_exp_val_cell(row_idx) = ap_exp_subj{il}(ie);
        end
      end
    end
  end

  % Trim preallocated arrays to actual size
  task_cell       = task_cell(1:row_idx);
  u_id_cell       = u_id_cell(1:row_idx);
  elec_cell       = elec_cell(1:row_idx);
  fooof_cell      = fooof_cell(1:row_idx);
  lat_cell        = lat_cell(1:row_idx);
  alpha_type_cell = alpha_type_cell(1:row_idx);
  gaze_meas_cell  = gaze_meas_cell(1:row_idx);
  bl_eeg_cell     = bl_eeg_cell(1:row_idx);
  bl_gaze_cell    = bl_gaze_cell(1:row_idx);
  subjectID_cell  = subjectID_cell(1:row_idx);
  Condition_cell  = Condition_cell(1:row_idx);
  alpha_val_cell  = alpha_val_cell(1:row_idx);
  gaze_val_cell   = gaze_val_cell(1:row_idx);
  n_trials_cell   = n_trials_cell(1:row_idx);
  ap_off_val_cell = ap_off_val_cell(1:row_idx);
  ap_exp_val_cell = ap_exp_val_cell(1:row_idx);

  tbl = table(task_cell, u_id_cell, elec_cell, fooof_cell, lat_cell, alpha_type_cell, gaze_meas_cell, ...
    bl_eeg_cell, bl_gaze_cell, subjectID_cell, Condition_cell, n_trials_cell, alpha_val_cell, gaze_val_cell, ...
    ap_off_val_cell, ap_exp_val_cell, ...
    'VariableNames', {'task', 'universe_id', 'electrodes', 'fooof', 'latency_ms', 'alpha_type', 'gaze_measure', ...
    'baseline_eeg', 'baseline_gaze', 'subjectID', 'Condition', 'n_trials_subject_condition', 'alpha', 'gaze_value', ...
    'aperiodic_offset', 'aperiodic_exponent'});

  if isempty(r2_sid)
    r2_tbl = table();
  else
    r2_tbl = table(r2_sid, r2_cond, r2_lat, r2_elec, ...
      r2_val, r2_err, r2_exp, r2_off, ...
      'VariableNames', {'subjectID', 'Condition', 'latency_ms', 'electrodes', ...
      'r_squared', 'fooof_error', 'aperiodic_exponent', 'aperiodic_offset'});
  end
end
