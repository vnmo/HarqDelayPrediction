%% ================= 공통 설정 =================
clear; clc;

% Monte Carlo 파라미터
iterLen  = 1e4;     % 필요시 1e5로 증가 가능
maxHarq  = 4;       % HARQ 1~4

% 시스템 파라미터 (MCS 14)
specEff      = 2.16;    % 5G LENA Table1 기준 MCS 14 spectral efficiency
cbs_target   = 816;     % 실제 사용 CBS
blockLen_fbl = cbs_target; % Theory(FBL)는 실제 CBS 사용
V            = 1.0;     % dispersion (네가 쓰는 값)
hvar         = 1.0;     % 평균 페이딩 전력

% SNR(dB) 축
snr_dB_vec = -10:1:20;
numSNR     = length(snr_dB_vec);

%% ================= 1. Theory (FBL-HARQ, CBS = 816) =================
bler_fbl = zeros(numSNR, maxHarq);

for k = 1:numSNR
    snr_lin = 10^(snr_dB_vec(k)/10);   % dB -> 선형
    
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

%% ================= 2. 5G LENA EESM 테이블 (MCS 14, CBS >= 816 중 가장 가까운 값) =================
% MCS 14의 CBS 후보들 (nr-eesm-t1.cc 기준)
cbs_lena_list = [3840, 4096, 4480, 4864, 5120, 5504, 6272, 6912, 7680, 8192];

% 816보다 크면서 가장 가까운 CBS 선택
idx      = find(cbs_lena_list >= cbs_target, 1, 'first');
cbs_lena = cbs_lena_list(idx);   % 이 경우 = 3840

% 선택된 CBS(3840)에 대한 LENA 테이블 값 (MCS 14, CBS 3840)
% nr-eesm-t1.cc 내용:
% SINR(dB): {8.0750, 8.6600, 8.9500, 9.2400, 9.8250, 10.7000}
% BLER    : {1, 0.9990385, 0.8827055, 0.3225703, 1e-4, 0}
snr_dB_lena = [8.0750, 8.6600, 8.9500, 9.2400, 9.8250, 10.7000];
bler_lena   = [1,      0.9990385, 0.8827055, 0.3225703, 1e-4, 0];

% log-도메인 보간 함수 (0,1 밖으로 안 나가게 클리핑)
bler_lena_fun = @(snr_dB) ...
    max(min( ...
        10.^interp1(snr_dB_lena, log10(bler_lena + 1e-12), ...
                    snr_dB, 'pchip', 'extrap'), ...
        1), 1e-12);

% LENA 기반 HARQ (EESM + Chase combining 가정: r회 전송 → 선형 SNR r배)
bler_lena_harq = zeros(numSNR, maxHarq);
for r = 1:maxHarq
    snr_eff_dB = snr_dB_vec + 10*log10(r);       % γ_eff,r(dB) = γ(dB) + 10log10(r)
    bler_lena_harq(:, r) = bler_lena_fun(snr_eff_dB);
end

%% ================= 3. 플롯: 색깔 매칭 + 선/점선 구분 =================
figure; hold on; grid on;

% HARQ 1~4용 색상: MATLAB 기본 'lines' 팔레트 사용
colors = lines(maxHarq);

for m = 1:maxHarq
    % Theory (FBL) – 실선
    semilogy(snr_dB_vec, bler_fbl(:, m), ...
        'LineStyle', '-', ...
        'Color', colors(m, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Theoretical Analysis HARQ %d Tx', m));
    
    % 5G LENA (EESM) – 점선, 동일 색상
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
