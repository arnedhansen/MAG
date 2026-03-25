%% AOC Multiverse — Trial-Level Data Preparation (Sternberg & N-Back)
% One-shot build: loads/computes EEG and gaze per trial for all
% electrode/FOOOF/latency/alpha/gaze/baseline combinations.
% Writes multiverse_sternberg.csv and multiverse_nback.csv.
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
% Model in R: alpha ~ gaze_value * Condition + (1|subjectID)

disp(upper('=== AOC MULTIVERSE PREP START ==='))

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
disp(upper(['Paths: base_data_in = ' base_data_in ', base_data_out = ' base_data_out ', base_features = ' base_features ', out_dir (CSV) = ' out_dir]))

if isempty(subjects)
    dirs = dir(base_features);
    dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {dirs.name};
    disp(upper(['Subject list from base_features: ' num2str(length(subjects)) ' folders.']))
end
if isempty(path_preproc)
    path_preproc = base_features;
    disp(upper('Using base_features as path_preproc (EEG and gaze under same root).'))
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

%% Get channel labels and indices (from first subject with power file)
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

%% Trial-level FOOOF estimator options
% Toggle between historical single-FFT and Welch-style trial spectra.
fooof_mode = 'welch';   % 'singleFFT' | 'welch'
fooof_mode_env = getenv('AOC_FOOOF_MODE');
if ~isempty(fooof_mode_env), fooof_mode = fooof_mode_env; end
welch_seg_len_sec = 0.5;      % 500 ms subsegments
welch_overlap = 0.5;          % 50% overlap
valid_fooof_modes = {'singleFFT', 'welch', 'BOTH'};
if ~ismember(fooof_mode, valid_fooof_modes)
  error('Invalid fooof_mode "%s". Valid options: %s', fooof_mode, strjoin(valid_fooof_modes, ', '));
end
% Parallel setup (Science Cloud safe): uses parfor when toolbox/pool is available.
use_parfor = init_parallel_pool_science_cloud();
if strcmpi(fooof_mode, 'BOTH')
  fooof_modes_to_run = {'singleFFT', 'welch'};
else
  fooof_modes_to_run = {fooof_mode};
end
write_untagged = numel(fooof_modes_to_run) == 1;
disp(upper(sprintf('Trial-level FOOOF mode request: %s | running: %s (seg=%.3fs, overlap=%.0f%%)', ...
  fooof_mode, strjoin(fooof_modes_to_run, ', '), welch_seg_len_sec, welch_overlap * 100)))
if ~write_untagged
  disp(upper('BOTH mode active: writing mode-tagged CSVs only (no untagged overwrite).'))
end

for im = 1:numel(fooof_modes_to_run)
  run_mode = fooof_modes_to_run{im};
  mode_tag = regexprep(run_mode, '[^A-Za-z0-9]+', '_');
  disp(upper(['=== RUN FOOOF MODE: ' run_mode ' ===']))

  %% Build Sternberg multiverse table
  disp(upper('--- STERNBERG TASK ---'))
  [tbl_s, r2_s] = build_task_multiverse('sternberg', subjects, path_preproc, base_features, ...
      ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
      baseline_eeg_opts, baseline_gaze_opts, ...
      n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange, ...
      run_mode, welch_seg_len_sec, welch_overlap, use_parfor);
  disp(upper(['Sternberg table rows: ' num2str(height(tbl_s))]))
  if height(tbl_s) == 0
    error('AOC_multiverse_prep:NoData', 'Sternberg table is empty.')
  end
  if write_untagged
    out_s_main = fullfile(out_dir, 'multiverse_sternberg.csv');
    writetable(tbl_s, out_s_main);
    disp(upper(['Written: ' out_s_main]))
  end
  out_s_mode = fullfile(out_dir, ['multiverse_sternberg_' mode_tag '.csv']);
  writetable(tbl_s, out_s_mode);
  disp(upper(['Written: ' out_s_mode]))
  if height(r2_s) > 0
    out_r2_s_mode = fullfile(r2_dir, ['fooof_r2_sternberg_' mode_tag '.csv']);
    if write_untagged
      out_r2_s_main = fullfile(r2_dir, 'fooof_r2_sternberg.csv');
      writetable(r2_s, out_r2_s_main);
      disp(upper(['Written: ' out_r2_s_main ' (' num2str(height(r2_s)) ' FOOOF fits)']))
    end
    writetable(r2_s, out_r2_s_mode);
    disp(upper(['Written: ' out_r2_s_mode]))
  end

  %% Build N-back multiverse table
  disp(upper('--- N-BACK TASK ---'))
  [tbl_n, r2_n] = build_task_multiverse('nback', subjects, path_preproc, base_features, ...
      ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
      baseline_eeg_opts, baseline_gaze_opts, ...
      n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange, ...
      run_mode, welch_seg_len_sec, welch_overlap, use_parfor);
  disp(upper(['N-back table rows: ' num2str(height(tbl_n))]))
  if height(tbl_n) == 0
    error('AOC_multiverse_prep:NoData', 'N-back table is empty.')
  end
  if write_untagged
    out_n_main = fullfile(out_dir, 'multiverse_nback.csv');
    writetable(tbl_n, out_n_main);
    disp(upper(['Written: ' out_n_main]))
  end
  out_n_mode = fullfile(out_dir, ['multiverse_nback_' mode_tag '.csv']);
  writetable(tbl_n, out_n_mode);
  disp(upper(['Written: ' out_n_mode]))
  if height(r2_n) > 0
    out_r2_n_mode = fullfile(r2_dir, ['fooof_r2_nback_' mode_tag '.csv']);
    if write_untagged
      out_r2_n_main = fullfile(r2_dir, 'fooof_r2_nback.csv');
      writetable(r2_n, out_r2_n_main);
      disp(upper(['Written: ' out_r2_n_main ' (' num2str(height(r2_n)) ' FOOOF fits)']))
    end
    writetable(r2_n, out_r2_n_mode);
    disp(upper(['Written: ' out_r2_n_mode]))
  end
end

disp(upper('=== AOC MULTIVERSE PREP DONE ==='))

%% ========== LOCAL FUNCTIONS ==========

function [idx_post, idx_occ] = get_channel_indices(labels)
  % Occipital: contains O or I; Parietal: contains P but NOT F; Posterior: union of occ + pari
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

function a = bandpower_trials(pow, chIdx, band)
  if isempty(chIdx), a = nan(size(pow.powspctrm, 1), 1); return, end
  f = pow.freq(:);
  bandIdx = f >= band(1) & f <= band(2);
  if sum(bandIdx) == 0, a = nan(size(pow.powspctrm, 1), 1); return, end
  x = pow.powspctrm(:, chIdx, bandIdx);
  a = squeeze(mean(mean(x, 2), 3));
  a = a(:);
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
  % Returns [gaze_dev, spl, vel, ms_cnt, bcea] at indices [1,2,3,4,5]
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
  % Baseline
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
  % Extract power from TFR for a given latency window: average over time → rpt_chan_freq
  cfg = []; cfg.latency = latency; cfg.trials = ind_cond;
  sel = ft_selectdata(cfg, tfr_all);
  sel.powspctrm = mean(sel.powspctrm, 4);
  if isfield(sel, 'time'), sel = rmfield(sel, 'time'); end
  sel.dimord = 'rpt_chan_freq';
  pow = sel;
end

function pow = compute_pow_cond_window(data_td, cond_inds, time_win, cfg_freq)
  % Compute trial-level power from time-domain EEG for one condition and time window
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

function pow_bl = compute_db_baseline_spectra(pow_task, pow_base)
  % Spectrum-level dB baseline: 10*log10(task / mean_baseline)
  pow_bl = pow_task;
  base_mean = mean(pow_base.powspctrm, 1);  % 1 x chan x freq
  % Align frequency axes if baseline has different resolution than task
  if length(pow_base.freq) ~= length(pow_task.freq) || max(abs(pow_base.freq - pow_task.freq)) > 1e-6
    base_aligned = nan(1, size(base_mean, 2), length(pow_task.freq));
    for ch = 1:size(base_mean, 2)
      base_aligned(1, ch, :) = interp1(pow_base.freq, squeeze(base_mean(1, ch, :)), pow_task.freq, 'linear', NaN);
    end
    base_mean = base_aligned;
  end
  base_mean(base_mean <= 0) = NaN;
  pow_bl.powspctrm = 10 * log10(bsxfun(@rdivide, pow_task.powspctrm, base_mean));
end

function pow_bl = compute_pct_baseline_spectra(pow_task, pow_base)
  % Spectrum-level percentage-change baseline: (task - base) / |base| * 100
  pow_bl = pow_task;
  base_mean = mean(pow_base.powspctrm, 1);  % 1 x chan x freq
  % Align frequency axes if baseline has different resolution than task
  if length(pow_base.freq) ~= length(pow_task.freq) || max(abs(pow_base.freq - pow_task.freq)) > 1e-6
    base_aligned = nan(1, size(base_mean, 2), length(pow_task.freq));
    for ch = 1:size(base_mean, 2)
      base_aligned(1, ch, :) = interp1(pow_base.freq, squeeze(base_mean(1, ch, :)), pow_task.freq, 'linear', NaN);
    end
    base_mean = base_aligned;
  end
  base_mean(base_mean <= 0) = NaN;
  pow_bl.powspctrm = bsxfun(@rdivide, ...
      bsxfun(@minus, pow_task.powspctrm, base_mean), abs(base_mean)) * 100;
end

function [f_osc, f_axis, ap_offset, ap_exponent, fooof_r2, fooof_err, fooof_n_segments] = ...
    run_fooof_from_raw(data_td, trial_idx, ch_labels, time_win, cfg_fooof, fooof_mode, welch_seg_len_sec, welch_overlap)
  % Run FOOOF from raw time-domain data: select trial, average ROI channels,
  % select time window, then call ft_freqanalysis_Arne_FOOOF.
  % Modes:
  %   - singleFFT: one periodogram from the full selected window
  %   - welch: split into 500 ms segments with 50% overlap and average before FOOOF
  % Also returns aperiodic parameters (offset, exponent) and fit quality (R², error).
  f_osc = []; f_axis = []; ap_offset = NaN; ap_exponent = NaN; fooof_r2 = NaN; fooof_err = NaN; fooof_n_segments = NaN;
  if isempty(data_td) || isempty(ch_labels), return, end
  try
    % Select single trial
    cfg_sel = []; cfg_sel.trials = trial_idx;
    d = ft_selectdata(cfg_sel, data_td);
    % Select time window
    cfg_lat = []; cfg_lat.latency = time_win;
    d = ft_selectdata(cfg_lat, d);
    % Select and average channels in ROI → single "ROI" channel
    cfg_ch = []; cfg_ch.channel = ch_labels; cfg_ch.avgoverchan = 'yes';
    d = ft_selectdata(cfg_ch, d);
    d.label = {'ROI'};
    % Build FOOOF input according to selected trial-level estimator mode
    if strcmp(fooof_mode, 'welch')
      [d_fooof, fooof_n_segments] = build_welch_segments_from_roi(d, welch_seg_len_sec, welch_overlap);
    else
      d_fooof = d;
      fooof_n_segments = 1;
    end
    % Run FOOOF (expects raw time-domain data, does FFT + FOOOF internally)
    if ~exist('ft_freqanalysis_Arne_FOOOF', 'file'), return, end
    out = ft_freqanalysis_Arne_FOOOF(cfg_fooof, d_fooof);
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
  catch
    f_osc = []; f_axis = [];
  end
end

function [d_seg, n_segments] = build_welch_segments_from_roi(d_roi, seg_len_sec, overlap_frac)
  % Convert one ROI trial into overlapped pseudo-trials so keeptrials='no'
  % yields a Welch-style averaged spectrum before FOOOF fitting.
  d_seg = d_roi;
  n_segments = 0;
  if isempty(d_roi) || ~isfield(d_roi, 'trial') || isempty(d_roi.trial) || isempty(d_roi.trial{1})
    return
  end
  x = double(d_roi.trial{1});
  if size(x,1) > 1
    x = mean(x, 1, 'omitnan');
  end
  if isfield(d_roi, 'fsample') && ~isempty(d_roi.fsample)
    fs = double(d_roi.fsample);
  else
    if ~isfield(d_roi, 'time') || isempty(d_roi.time) || numel(d_roi.time{1}) < 2
      return
    end
    fs = 1 / median(diff(double(d_roi.time{1})));
  end
  seg_n = max(2, round(seg_len_sec * fs));
  hop_n = max(1, round(seg_n * (1 - overlap_frac)));
  n = numel(x);
  if n < seg_n
    d_seg = d_roi;
    n_segments = 1;
    return
  end
  starts = 1:hop_n:(n - seg_n + 1);
  n_segments = numel(starts);
  d_seg.trial = cell(1, n_segments);
  d_seg.time = cell(1, n_segments);
  d_seg.sampleinfo = nan(n_segments, 2);
  for k = 1:n_segments
    idx = starts(k):(starts(k) + seg_n - 1);
    d_seg.trial{k} = x(idx);
    if isfield(d_roi, 'time') && ~isempty(d_roi.time) && ~isempty(d_roi.time{1})
      t0 = d_roi.time{1}(idx(1));
      d_seg.time{k} = t0 + (0:seg_n-1) ./ fs;
    else
      d_seg.time{k} = (0:seg_n-1) ./ fs;
    end
    d_seg.sampleinfo(k,:) = [1 seg_n];
  end
  d_seg.label = {'ROI'};
  d_seg.fsample = fs;
  if isfield(d_seg, 'trialinfo')
    if isempty(d_seg.trialinfo)
      d_seg = rmfield(d_seg, 'trialinfo');
    else
      d_seg.trialinfo = repmat(d_seg.trialinfo(1, :), n_segments, 1);
    end
  end
end

function use_parfor = init_parallel_pool_science_cloud()
  use_parfor = false;
  try
    has_parallel = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
  catch
    has_parallel = false;
  end
  if ~has_parallel
    disp(upper('PARFOR DISABLED: PARALLEL TOOLBOX NOT AVAILABLE.'))
    return
  end
  try
    p = gcp('nocreate');
    if isempty(p)
      parpool('local');
      disp(upper('PARFOR ENABLED: STARTED LOCAL PARPOOL.'))
    else
      disp(upper(['PARFOR ENABLED: USING EXISTING POOL (' num2str(p.NumWorkers) ' WORKERS).']))
    end
    use_parfor = true;
  catch ME
    disp(upper(['PARFOR DISABLED: ' ME.message]))
    use_parfor = false;
  end
end

%% ========== MAIN BUILD FUNCTION ==========

function [tbl, r2_tbl] = build_task_multiverse(task_name, subjects, path_preproc, base_features, ...
    ch_sets, ch_label_sets, electrodes_opts, fooof_opts, latency_opts, alpha_opts, gaze_opts, gaze_col_map, ...
    baseline_eeg_opts, baseline_gaze_opts, ...
    n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze, n_universes, alphaRange, ...
    fooof_mode, welch_seg_len_sec, welch_overlap, use_parfor)

  disp(upper(['Building multiverse table for task: ' task_name ' (' num2str(n_universes) ' universes)']))
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

  % Freq analysis config (Hanning)
  cfg_hann = []; cfg_hann.method = 'mtmfft'; cfg_hann.taper = 'hanning';
  cfg_hann.foilim = [2 40]; cfg_hann.pad = 5; cfg_hann.keeptrials = 'yes';

  % FOOOF config (Hanning)
  cfg_fooof = []; cfg_fooof.method = 'mtmfft'; cfg_fooof.taper = 'hanning';
  cfg_fooof.foilim = [2 40]; cfg_fooof.pad = 5; cfg_fooof.output = 'fooof'; cfg_fooof.keeptrials = 'no';

  % Latency windows
  lat_windows = {[0 0.5], [0 1], [0 2], [1 2]};
  lat_durs    = [0.5, 1, 2, 1];
  t_base = [-0.5 -0.25];
  fsample = 500;

  % Preallocate output arrays
  task_cell = {}; u_id_cell = []; elec_cell = {}; fooof_cell = {}; lat_cell = {};
  alpha_type_cell = {}; gaze_meas_cell = {}; bl_eeg_cell = {}; bl_gaze_cell = {};
  subjectID_cell = []; Trial_cell = []; Condition_cell = [];
  alpha_val_cell = []; gaze_val_cell = [];
  ap_off_val_cell = []; ap_exp_val_cell = [];

  % FOOOF R² collection (compact: one row per FOOOF call)
  r2_sid = []; r2_cond = []; r2_trial = []; r2_lat = {}; r2_elec = {};
  r2_mode = {}; r2_nseg = [];
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

    %% ====== EEG: Load pre-computed power (Hanning, raw) ======
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
    % convolution window causes edge NaN at t≥1.75s. mean() over time then
    % propagates NaN to the entire trial-averaged spectrum. If all NaN,
    % clear so the fallback (direct FFT on time-domain data) kicks in.
    for c = 1:length(pow_full)
      if ~isempty(pow_full{c}) && isfield(pow_full{c}, 'powspctrm') && ...
          all(isnan(pow_full{c}.powspctrm(:)))
        disp(upper(['  WARNING: pow_full{' num2str(c) '} is all NaN (TFR edge effect). Will recompute from time-domain.']))
        pow_full{c} = [];
      end
    end

    %% ====== Load TFR (Hanning, for arbitrary windows) ======
    tfr_path = fullfile(eeg_dir, tfr_file);
    eeg_tfr_path = fullfile(eeg_dir, eeg_tfr_file);
    eeg_td_path = fullfile(eeg_dir, eeg_td_file);

    tfr_loaded = [];
    if isfile(tfr_path)
      disp(upper('  Loading precomputed TFR (Hanning).'))
      load(tfr_path, 'tfr_all');
      tfr_loaded = tfr_all;
      disp(upper('  TFR loaded.'))
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
      disp(upper('  WARNING: No TFR or time-domain EEG. Hanning 0-500ms/1-2s and all FOOOFed will be NaN.'))
    end
    if isempty(data_td)
      disp(upper('  WARNING: No time-domain EEG loaded → all FOOOFed alpha universes will be NaN for this subject.'))
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
      trials_list = trial_num(trl_idx);

      % EEG-ET trial alignment
      n_eeg_trials = size(pow_early{cond}.powspctrm, 1);
      if n_trials_cond ~= n_eeg_trials
        disp(upper(['  WARNING: ET=' num2str(n_trials_cond) ' vs EEG=' num2str(n_eeg_trials) ' for cond ' num2str(cval) '. Using min.']))
        n_trials_cond = min(n_trials_cond, n_eeg_trials);
        trl_idx = trl_idx(1:n_trials_cond);
        trials_list = trial_num(trl_idx);
      end
      disp(upper(['  Condition ' num2str(cval) ': ' num2str(n_trials_cond) ' trials.']))

      %% ====== Build power struct grid: pow_s{lat}, pow_db_s{lat}, pow_pct_s{lat} ======
      pow_s     = cell(4, 1);  % raw
      pow_db_s  = cell(4, 1);  % dB-baselined
      pow_pct_s = cell(4, 1);  % pct-change-baselined

      % --- Pre-computed Hanning (0-1s, 0-2s) ---
      pow_s{2} = pow_early{cond};       % 0-1s, raw
      pow_s{3} = pow_full{cond};        % 0-2s, raw

      % --- TFR-based Hanning (0-500ms, 1-2s) ---
      if ~isempty(tfr_loaded)
        cond_inds_tfr = find(tfr_loaded.trialinfo(:,1) == cond_codes(cond));
        if isempty(pow_s{1})
          pow_s{1} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{1});
        end
        if isempty(pow_s{4})
          pow_s{4} = get_power_from_TFR(tfr_loaded, cond_inds_tfr, lat_windows{4});
        end
      end

      % --- Time-domain fill (remaining Hanning gaps) ---
      cond_inds_td = [];  % initialize before block; needed for FOOOF later
      if ~isempty(data_td)
        cond_inds_td = find(data_td.trialinfo(:,1) == cond_codes(cond));
        for il = 1:4
          if isempty(pow_s{il})
            disp(upper(['    Computing power from time-domain: lat=' num2str(il)]))
            pow_s{il} = compute_pow_cond_window(data_td, cond_inds_td, lat_windows{il}, cfg_hann);
          end
        end
      end

      % --- Baseline spectrum for dB and percentage-change correction ---
      pow_base = [];
      if ~isempty(tfr_loaded)
        try
          pow_base = get_power_from_TFR(tfr_loaded, cond_inds_tfr, t_base);
        catch
          disp(upper('    WARNING: TFR baseline extraction failed (baseline window outside TFR range). Trying time-domain.'))
          pow_base = [];
        end
      end
      if isempty(pow_base) && ~isempty(data_td) && ~isempty(cond_inds_td)
        pow_base = compute_pow_cond_window(data_td, cond_inds_td, t_base, cfg_hann);
      end
      for il = 1:4
        if ~isempty(pow_s{il}) && ~isempty(pow_base)
          if length(pow_s{il}.freq) ~= length(pow_base.freq)
            disp(upper(['    NOTE: Freq axis mismatch for lat=' num2str(il) ...
              ' (task=' num2str(length(pow_s{il}.freq)) ' vs base=' num2str(length(pow_base.freq)) ...
              ' bins). Interpolating baseline.']))
          end
          pow_db_s{il}  = compute_db_baseline_spectra(pow_s{il}, pow_base);
          pow_pct_s{il} = compute_pct_baseline_spectra(pow_s{il}, pow_base);
        end
      end

      %% ====== Alpha: cell array {lat, bl, fo}(trial, elec_alpha_col) ======
      % fo: 1=FOOOFed, 2=nonFOOOFed
      n_cols = n_elec * n_alpha;
      alpha = cell(n_lat, n_bl_eeg, n_fooof);
      for il = 1:n_lat; for ibl = 1:n_bl_eeg; for ifo = 1:n_fooof
        alpha{il, ibl, ifo} = nan(n_trials_cond, n_cols);
      end; end; end

      % Aperiodic parameters from FOOOF: {lat}(trial, elec)
      ap_offset_arr = cell(n_lat, 1);
      ap_exp_arr = cell(n_lat, 1);
      for il = 1:n_lat
        ap_offset_arr{il} = nan(n_trials_cond, n_elec);
        ap_exp_arr{il} = nan(n_trials_cond, n_elec);
      end

      disp(upper('  Computing alpha (raw + FOOOF × raw/dB/pct_change) per trial.'))
      for il = 1:n_lat
        % --- nonFOOOFed: bandpower from pre-computed / computed power structs ---
        for ibl = 1:n_bl_eeg
          if ibl == 1, p = pow_s{il}; elseif ibl == 2, p = pow_db_s{il}; else, p = pow_pct_s{il}; end
          if isempty(p) || ~isfield(p, 'powspctrm'), continue; end
          for ie = 1:n_elec
            chIdx = ch_sets{ie};
            for ia = 1:n_alpha
              col = (ie-1)*n_alpha + ia;
              if ia == 1, band = alphaRange; else, band = IAF_band; end
              alpha{il, ibl, 2}(:, col) = bandpower_trials(p, chIdx, band);
            end
          end
        end

        % --- FOOOFed: from raw time-domain data ---
        % ft_freqanalysis_Arne_FOOOF needs raw (time-domain) input, NOT a freq struct.
        % It internally performs FFT + FOOOF. We pass single-trial, ROI-averaged data.
        % FOOOF result is baseline-independent → same value for raw, dB, and pct_change baseline.
        if ~isempty(data_td) && ~isempty(cond_inds_td)
          for ie = 1:n_elec
            ch_lbl = ch_label_sets{ie};
            n_td = min(n_trials_cond, length(cond_inds_td));
            ap_off_vec = nan(n_td, 1); ap_exp_vec = nan(n_td, 1);
            f_r2_vec = nan(n_td, 1); f_err_vec = nan(n_td, 1); f_nseg_vec = nan(n_td, 1);
            fooof_alpha_vec = nan(n_td, n_alpha);
            if use_parfor
              parfor tr = 1:n_td
                td_idx = cond_inds_td(tr);
                [f_osc, f_axis, ap_off, ap_exp, f_r2, f_err, f_nseg] = run_fooof_from_raw( ...
                  data_td, td_idx, ch_lbl, lat_windows{il}, cfg_fooof, fooof_mode, welch_seg_len_sec, welch_overlap);
                ap_off_vec(tr) = ap_off; ap_exp_vec(tr) = ap_exp;
                f_r2_vec(tr) = f_r2; f_err_vec(tr) = f_err; f_nseg_vec(tr) = f_nseg;
                if isempty(f_osc), continue, end
                for ia = 1:n_alpha
                  if ia == 1, band = alphaRange; else, band = IAF_band; end
                  bandIdx = f_axis >= band(1) & f_axis <= band(2);
                  if any(bandIdx)
                    fooof_alpha_vec(tr, ia) = mean(f_osc(bandIdx), 'omitnan');
                  end
                end
              end
            else
              for tr = 1:n_td
                td_idx = cond_inds_td(tr);
                [f_osc, f_axis, ap_off, ap_exp, f_r2, f_err, f_nseg] = run_fooof_from_raw( ...
                  data_td, td_idx, ch_lbl, lat_windows{il}, cfg_fooof, fooof_mode, welch_seg_len_sec, welch_overlap);
                ap_off_vec(tr) = ap_off; ap_exp_vec(tr) = ap_exp;
                f_r2_vec(tr) = f_r2; f_err_vec(tr) = f_err; f_nseg_vec(tr) = f_nseg;
                if isempty(f_osc), continue, end
                for ia = 1:n_alpha
                  if ia == 1, band = alphaRange; else, band = IAF_band; end
                  bandIdx = f_axis >= band(1) & f_axis <= band(2);
                  if any(bandIdx)
                    fooof_alpha_vec(tr, ia) = mean(f_osc(bandIdx), 'omitnan');
                  end
                end
              end
            end
            for tr = 1:n_td
              ap_offset_arr{il}(tr, ie) = ap_off_vec(tr);
              ap_exp_arr{il}(tr, ie) = ap_exp_vec(tr);
              for ia = 1:n_alpha
                col = (ie-1)*n_alpha + ia;
                if isfinite(fooof_alpha_vec(tr, ia))
                  for ibl = 1:n_bl_eeg
                    alpha{il, ibl, 1}(tr, col) = fooof_alpha_vec(tr, ia);
                  end
                end
              end
              if isfinite(f_r2_vec(tr)) || isfinite(f_err_vec(tr)) || isfinite(ap_exp_vec(tr)) || isfinite(ap_off_vec(tr))
                r2_sid(end+1,1) = sid_num; r2_cond(end+1,1) = cval;
                r2_trial(end+1,1) = trials_list(tr);
                r2_lat{end+1,1} = latency_opts{il}; r2_elec{end+1,1} = electrodes_opts{ie};
                r2_mode{end+1,1} = fooof_mode; r2_nseg(end+1,1) = f_nseg_vec(tr);
                r2_val(end+1,1) = f_r2_vec(tr); r2_err(end+1,1) = f_err_vec(tr);
                r2_exp(end+1,1) = ap_exp_vec(tr); r2_off(end+1,1) = ap_off_vec(tr);
              end
            end
          end
        else
          disp(upper('    WARNING: No time-domain EEG data → FOOOFed alpha will be NaN for this condition.'))
        end
      end

      %% ====== Gaze: 4 metrics × 4 windows × raw + dB + pct_change ======
      disp(upper('  Computing gaze (dev, SPL, vel, MS, BCEA) x 4 windows x raw + dB + pct_change.'))
      gaze_raw_all = cell(n_lat, 1); gaze_db_all = cell(n_lat, 1); gaze_pct_all = cell(n_lat, 1);
      for il = 1:n_lat
        gaze_raw_all{il} = nan(n_trials_cond, 5);
        gaze_db_all{il}  = nan(n_trials_cond, 5);
        gaze_pct_all{il} = nan(n_trials_cond, 5);
      end
      for tr = 1:n_trials_cond
        trl_glob = trl_idx(tr);
        raw_et = dataET.trial{trl_glob};
        t = dataET.time{trl_glob};
        if size(raw_et,1) < 2, continue; end
        x = raw_et(1,:); y = raw_et(2,:);
        if length(t) ~= length(x), continue; end
        for il = 1:n_lat
          [gaze_raw_all{il}(tr,:), gaze_db_all{il}(tr,:), gaze_pct_all{il}(tr,:)] = ...
            compute_gaze_one_window(x, y, t, lat_windows{il}, lat_durs(il), fsample, t_base);
        end
      end

      %% ====== Append rows: trial × universe ======
      disp(upper(['  Appending ' num2str(n_trials_cond * n_universes) ' rows for subject ' sid ' cond ' num2str(cval) '.']))
      for tr = 1:n_trials_cond
        for u = 1:n_universes
          [ie, ifo, il, ia, ig, ibe, ibg] = ind2sub( ...
            [n_elec, n_fooof, n_lat, n_alpha, n_gaze, n_bl_eeg, n_bl_gaze], u);
          task_cell{end+1, 1}       = task_name;
          u_id_cell(end+1, 1)       = u;
          elec_cell{end+1, 1}       = electrodes_opts{ie};
          fooof_cell{end+1, 1}      = fooof_opts{ifo};
          lat_cell{end+1, 1}        = latency_opts{il};
          alpha_type_cell{end+1, 1} = alpha_opts{ia};
          gaze_meas_cell{end+1, 1}  = gaze_opts{ig};
          bl_eeg_cell{end+1, 1}     = baseline_eeg_opts{ibe};
          bl_gaze_cell{end+1, 1}    = baseline_gaze_opts{ibg};
          subjectID_cell(end+1, 1)  = sid_num;
          Trial_cell(end+1, 1)      = trials_list(tr);
          Condition_cell(end+1, 1)  = cval;

          % Alpha lookup
          col_alpha = (ie-1)*n_alpha + ia;
          av = alpha{il, ibe, ifo}(tr, col_alpha);
          alpha_val_cell(end+1, 1) = av;

          % Gaze lookup (gaze_col_map maps gaze_opts index to column in 5-element gaze array)
          gcol = gaze_col_map(ig);
          if ibg == 1
            gv = gaze_raw_all{il}(tr, gcol);
          elseif ibg == 2
            gv = gaze_db_all{il}(tr, gcol);
          else
            gv = gaze_pct_all{il}(tr, gcol);
          end
          gaze_val_cell(end+1, 1) = gv;

          % Aperiodic lookup (only valid for FOOOFed universes)
          if ifo == 1
            ap_off_val_cell(end+1, 1) = ap_offset_arr{il}(tr, ie);
            ap_exp_val_cell(end+1, 1) = ap_exp_arr{il}(tr, ie);
          else
            ap_off_val_cell(end+1, 1) = NaN;
            ap_exp_val_cell(end+1, 1) = NaN;
          end
        end
      end
    end
  end

  tbl = table(task_cell, u_id_cell, elec_cell, fooof_cell, lat_cell, alpha_type_cell, gaze_meas_cell, ...
    bl_eeg_cell, bl_gaze_cell, subjectID_cell, Trial_cell, Condition_cell, alpha_val_cell, gaze_val_cell, ...
    ap_off_val_cell, ap_exp_val_cell, ...
    'VariableNames', {'task', 'universe_id', 'electrodes', 'fooof', 'latency_ms', 'alpha_type', 'gaze_measure', ...
    'baseline_eeg', 'baseline_gaze', 'subjectID', 'Trial', 'Condition', 'alpha', 'gaze_value', ...
    'aperiodic_offset', 'aperiodic_exponent'});

  if isempty(r2_sid)
    r2_tbl = table();
  else
    r2_tbl = table(r2_sid, r2_cond, r2_trial, r2_lat, r2_elec, r2_mode, r2_nseg, ...
      r2_val, r2_err, r2_exp, r2_off, ...
      'VariableNames', {'subjectID', 'Condition', 'Trial', 'latency_ms', 'electrodes', ...
      'fooof_mode', 'welch_n_segments', 'r_squared', 'fooof_error', 'aperiodic_exponent', 'aperiodic_offset'});
  end
end
