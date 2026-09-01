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

targetMcsList = [20 14 12 9 5];

if ~exist("testbedresults", "dir")
    mkdir("testbedresults");
end

figure(1);

for dIter = 1:length(delayTargetArray)
    d = delayTargetArray(dIter);

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

    NrbPlot = nan(1, length(NrbList));
    McsPlot = nan(1, length(NrbList));

    McsPrev = -1;

    for nrbIter = 1:length(NrbList)
        Nrb = NrbList(nrbIter);
        Nre = Nrb * 134.16;
        B = 134e3 * Nrb; %#ok<NASGU>

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

        if ismember(mcs, targetMcsList) && (mcs ~= McsPrev)

            [~, ~, eta_target] = parseMcs(mcs, 1);
            Nrb_target = ceil(pktLen / (134.16 * eta_target));

            NrbPlot(nrbIter) = Nrb_target;
            McsPlot(nrbIter) = mcs;

            p_Harq = getLenaHarqBler( ...
                simLen_MonteCarlo, ...
                mcs, ...
                SNRdb_avgRcvd, ...
                hvar, ...
                M_Harq, ...
                2, ...
                816, ...
                Nrb_target);

            expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
            for i = 1:M_Harq-1
                expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) ...
                    + i * prod(p_Harq(1:i-1)) * (1 - p_Harq(i));
            end

            dvp_Harq(nrbIter) = getDvpHarqTestbed( ...
                d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, mcs);

            McsPrev = mcs;
        end
    end

    dvp_Harq(dvp_Harq > 1) = 1;

    if isSimRequired
    else
        validIdx = ~isnan(dvp_Harq) & ~isnan(NrbPlot) & ~isnan(McsPlot);

        NrbList_valid = NrbPlot(validIdx);
        dvp_Harq_valid = dvp_Harq(validIdx);
        Mcs_valid = McsPlot(validIdx);

        [NrbList_valid, sortIdx] = sort(NrbList_valid);
        dvp_Harq_valid = dvp_Harq_valid(sortIdx);
        Mcs_valid = Mcs_valid(sortIdx);

        writematrix(NrbList_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "NRBListTheory.csv");
        writematrix(dvp_Harq_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "DVPListTheory.csv");
        writematrix(Mcs_valid, "testbedresults/" + string(delayTargetArrayOrig(dIter)) + "MCSListTheory.csv");

        plot(NrbList_valid, log10(dvp_Harq_valid), lineOrder(1), ...
            'Marker', markerOrder(dIter), ...
            'MarkerSize', 10, ...
            'color', corder(dIter,:), ...
            'LineWidth', 3);
        hold on;

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

