%% ================= Common settings =================
clear; clc;

% Monte Carlo parameters
iterLen  = 1e4;     % Increase to 1e5 if needed
maxHarq  = 4;       % HARQ 1~4

% System parameters (MCS 14)
specEff      = 2.16;    % MCS 14 spectral efficiency from 5G LENA Table 1
cbs_target   = 816;     % CBS used in this evaluation
blockLen_fbl = cbs_target; % Use the configured CBS for FBL theory
V            = 1.0;     % Configured channel dispersion
hvar         = 1.0;     % Mean fading power

% SNR axis in dB
snr_dB_vec = -10:1:20;
numSNR     = length(snr_dB_vec);

%% ================= 1. Theory (FBL-HARQ, CBS = 816) =================
bler_fbl = zeros(numSNR, maxHarq);

for k = 1:numSNR
    snr_lin = 10^(snr_dB_vec(k)/10);   % dB -> linear
    
    expectedBlerVector = getFblPerVector_Harq( ...
        iterLen, ...
        snr_lin, ...
        specEff, ...
        blockLen_fbl, ...
        V, ...
        hvar, ...
        maxHarq ...
    );
    
    bler_fbl(k, :) = expectedBlerVector;   % 1 x maxHarq
end

%% ================= 2. 5G LENA EESM table (MCS 14, smallest CBS >= 816) =================
% Available CBS values for MCS 14 (from nr-eesm-t1.cc)
cbs_lena_list = [3840, 4096, 4480, 4864, 5120, 5504, 6272, 6912, 7680, 8192];

% Select the smallest available CBS at least as large as the target
idx      = find(cbs_lena_list >= cbs_target, 1, 'first');
cbs_lena = cbs_lena_list(idx);   % 3840 in this case

% LENA table values for the selected CBS (MCS 14, CBS 3840)
% Values from nr-eesm-t1.cc:
% SINR(dB): {8.0750, 8.6600, 8.9500, 9.2400, 9.8250, 10.7000}
% BLER    : {1, 0.9990385, 0.8827055, 0.3225703, 1e-4, 0}
snr_dB_lena = [8.0750, 8.6600, 8.9500, 9.2400, 9.8250, 10.7000];
bler_lena   = [1,      0.9990385, 0.8827055, 0.3225703, 1e-4, 0];

% Log-domain interpolation, clipped to the probability range
bler_lena_fun = @(snr_dB) ...
    max(min( ...
        10.^interp1(snr_dB_lena, log10(bler_lena + 1e-12), ...
                    snr_dB, 'pchip', 'extrap'), ...
        1), 1e-12);

% LENA-based HARQ (EESM + Chase combining: r attempts scale linear SNR by r)
bler_lena_harq = zeros(numSNR, maxHarq);
for r = 1:maxHarq
    snr_eff_dB = snr_dB_vec + 10*log10(r);       % γ_eff,r(dB) = γ(dB) + 10log10(r)
    bler_lena_harq(:, r) = bler_lena_fun(snr_eff_dB);
end

%% ================= 3. Plot: matched colors with solid/dashed lines =================
figure; hold on; grid on;

% Use the MATLAB 'lines' palette for HARQ attempts 1:4
colors = lines(maxHarq);

for m = 1:maxHarq
    % Theory (FBL) - solid line
    semilogy(snr_dB_vec, bler_fbl(:, m), ...
        'LineStyle', '-', ...
        'Color', colors(m, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Theoretical Analysis HARQ %d Tx', m));
    
    % 5G LENA (EESM) - dashed line, same color
    semilogy(snr_dB_vec, bler_lena_harq(:, m), ...
        'LineStyle', '--', ...
        'Color', colors(m, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('5G LENA HARQ %d Tx', m));
end

xlabel('SNR (dB)');
ylabel('BLER');
title(sprintf('MCS 14: Theory  vs 5G LENA '));
legend('Location', 'southwest');
ylim([1e-4 1]);
xlim([snr_dB_vec(1) snr_dB_vec(end)]);
grid on;
set(gca, 'YMinorGrid', 'on');
