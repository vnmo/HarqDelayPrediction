close all;
clearvars -except h1 h2 h3; clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

delayTargetArrayOrig = [2 3 5 6 8 10];

isReTxPrioritised = true;
isSimRequired = false;
simLen = 1e6;
simLen_MonteCarlo = 1e6;

Qmax = Inf;
M = Inf;
Qmax_Harq = 16;
M_Harq = 4;

T = 1e-3;
hvar = 1; % mean of h^2
V = 1;

SNRdb_avgRcvd = 10;

snr_avgRcvd = 10.^(SNRdb_avgRcvd/10); %#ok<NASGU>
f_inverse = 3;
f = 1/f_inverse;

delta_Orig = 1; % Feedback delay
pktLen = 8 * 108; % Packet length in bits
CSO = 0; % Control Signalling Overhead

zeta = 1; % Decoding delay
delta = delta_Orig + zeta;
delayTargetArray = delayTargetArrayOrig * 1e-3;

if min(delayTargetArray) < 1e-3
    error("Invalid delayTarget,zeta combination. Correct it!")
end

% ===== 사용할 MCS만 지정 =====
targetMcsList = [20 14 12 9 5];

% 출력 폴더
if ~exist("testbedresults", "dir")
    mkdir("testbedresults");
end

figure(1);

for dIter = 1:length(delayTargetArray)
    d = delayTargetArray(dIter);

    % eta is spectral efficiency
    [~, ~, etaMax] = parseMcs(28,1);
    Nrb_min = ceil(pktLen / (134.16 * etaMax));
    [~, ~, etaMin] = parseMcs(0,1);
    Nrb_max = ceil(pktLen / (134.16 * etaMin));

    NrbList = Nrb_min:1:Nrb_max;

    expNumTx_Arq = zeros(1, length(NrbList)); %#ok<NASGU>
    dvp_Arq = nan(1, length(NrbList)); %#ok<NASGU>
    dvp_Arq_sim = nan(1, length(NrbList)); %#ok<NASGU>
    expNumTx_Harq = zeros(1, length(NrbList));
    dvp_Harq = nan(1, length(NrbList));
    dvp_Harq_sim = nan(1, length(NrbList)); %#ok<NASGU>

    % 실제 plot/save에 쓸 대표 NRB와 MCS
    NrbPlot = nan(1, length(NrbList));
    McsPlot = nan(1, length(NrbList));

    % 같은 MCS를 여러 번 저장하지 않기 위한 변수
    McsPrev = -1;

    for nrbIter = 1:length(NrbList)
        Nrb = NrbList(nrbIter);
        Nre = Nrb * 134.16;
        B = 134e3 * Nrb; %#ok<NASGU>

        % 현재 Nrb에서 packet을 보낼 수 있는 최소 MCS 선택
        etaMin_local = pktLen / (134.16 * Nrb);
        eta = [];
        mcs = [];

        for mcsItem = 1:28
            [~, ~, eta_mcsItem] = parseMcs(mcsItem, 1);
            if eta_mcsItem >= etaMin_local
                eta = eta_mcsItem;
                mcs = mcsItem;
                break;
            end
        end

        disp(mcs);
        dataBits = Nre * eta; %#ok<NASGU>

        % ===== target MCS만 사용, 같은 MCS는 첫 등장 1번만 계산 =====
        if ismember(mcs, targetMcsList) && (mcs ~= McsPrev)

            % x축용 대표 PRB per slot은 현재 scan Nrb가 아니라
            % 해당 MCS의 spectral efficiency로 직접 계산
            [~, ~, eta_target] = parseMcs(mcs, 1);
            Nrb_target = ceil(pktLen / (134.16 * eta_target));

            NrbPlot(nrbIter) = Nrb_target;
            McsPlot(nrbIter) = mcs;

            %% ========= HARQ ================================================
            % p_Harq = getFblPerVector_Harq(simLen_MonteCarlo, snr_avgRcvd, eta_target, Nrb_target*134.16, V, hvar, M_Harq);
            p_Harq = getLenaHarqBler( ...
                simLen_MonteCarlo, ...
                mcs, ...
                SNRdb_avgRcvd, ...
                hvar, ...
                M_Harq, ...
                2, ...
                816, ...
                Nrb_target);

            %% Expected #transmissions
            expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
            for i = 1:M_Harq-1
                expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) ...
                    + i * prod(p_Harq(1:i-1)) * (1 - p_Harq(i));
            end

            %% DVP with HARQ using testbed queue PMF
            % getDvpHarqTestbed 내부에서 comparison_queue_pmf.csv의
            % resMCS{mcs} column을 읽는다고 가정
            dvp_Harq(nrbIter) = getDvpHarqTestbed( ...
                d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, mcs);

            % 같은 MCS 반복 저장 방지
            McsPrev = mcs;
        end
    end

    dvp_Harq(dvp_Harq > 1) = 1;

    if isSimRequired
        % 필요시 시뮬레이션 plot 추가
    else
        % ===== 선택된 점만 추출 =====
        validIdx = ~isnan(dvp_Harq) & ~isnan(NrbPlot) & ~isnan(McsPlot);

        NrbList_valid = NrbPlot(validIdx);
        dvp_Harq_valid = dvp_Harq(validIdx);
        Mcs_valid = McsPlot(validIdx);

        % x축 순서대로 정렬
        [NrbList_valid, sortIdx] = sort(NrbList_valid);
        dvp_Harq_valid = dvp_Harq_valid(sortIdx);
        Mcs_valid = Mcs_valid(sortIdx);

        % 저장
        writematrix(NrbList_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "NRBListTheory.csv");
        writematrix(dvp_Harq_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "DVPListTheory.csv");
        writematrix(Mcs_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "MCSListTheory.csv");

        % plot
        plot(NrbList_valid, log10(dvp_Harq_valid), lineOrder(1), ...
            'Marker', markerOrder(dIter), ...
            'MarkerSize', 10, ...
            'color', corder(dIter,:), ...
            'LineWidth', 3);
        hold on;

        % 각 점에 MCS 라벨 표시하고 싶으면 이 블록 주석 해제
        % for k = 1:length(NrbList_valid)
        %     text(NrbList_valid(k) + 0.15, log10(dvp_Harq_valid(k)), ...
        %         "MCS " + string(Mcs_valid(k)), ...
        %         'FontSize', 12, 'Color', corder(dIter,:));
        % end
    end
end

figure(1)
set(gca, 'FontSize', 28)
grid on;

legend( ...
    'HARQ: $d = 2$ ms', ...
    'HARQ: $d = 3$ ms', ...
    'HARQ: $d = 5$ ms', ...
    'HARQ: $d = 6$ ms', ...
    'HARQ: $d = 8$ ms', ...
    'HARQ: $d = 10$ ms', ...
    'fontsize', 32, ...
    'interpreter', 'latex', ...
    'location', 'southwest', ...
    'NumColumns', 2, ...
    'Orientation', 'horizontal');

xlabel('PRB per slot', 'fontsize', 42, 'interpreter', 'latex');
ylabel('DVP', 'fontsize', 42, 'interpreter', 'latex');

ylim([-6, 0])
xlim([0, 20])

yticks([-6, -5, -4, -3, -2, -1, 0])
yticklabels({'10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '1'})

% saveas(gcf, 'Figures/blerVsNrb_vsDelayTarget.fig')
% saveas(gcf, 'Figures/blerVsNrb_vsDelayTarget.epsc')
% movefile('Figures/blerVsNrb_vsDelayTarget.epsc', 'Figures/blerVsNrb_vsDelayTarget.eps')