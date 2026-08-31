% Queueing simulation: DVP sweeps using simulate_queue_IR_cso.
% Outputs are *Simulation.csv; getDvpHarqCSO/getQLenPmfCSO are not called.
% The simulation uses packet lists and a slot clock for periodic CSO.
% k_xi is the CSP length; phase is mod(slotId,k_xi) inside the simulator.
% Scenario settings are defined below. simLen is a slot count.
close all;
clearvars; clc;

%% Sweep settings
delayTargetArrayOrig = [2 8];   % in ms
f_inverse_list = [3];%[3 6 8 10];
sinr_list = [20];%[10 20 ]
k_xi_list = [80];%[80 320]

%% Simulation options
isReTxPrioritised = true;
simLen = 1e7;
maxQueueSize = 16;
isDetArr = false;

%% HARQ settings
M_Harq = 4;
rttSlots = 3;                    % 1 Tx + 1 decoding + feedback
feedbackDelay = rttSlots - 2;    % slots

%% System settings
T = 1e-3;            % slot duration in seconds
pktLen = 8 * 108;    % packet length in bits
hvar = 1;

for sinr_index = 1:length(sinr_list)
    for par_index = 1:length(f_inverse_list)
        for kxi_index = 1:length(k_xi_list)

            %% Current sweep point
            SNRdb_avgRcvd = sinr_list(sinr_index);
            f_inverse = f_inverse_list(par_index);
            f = 1 / f_inverse;
            k_xi = k_xi_list(kxi_index);

            saving_folder = "sinr" + string(SNRdb_avgRcvd) + ...
                            "_par" + string(f_inverse) + ...
                            "_srs" + string(k_xi);
            mkdir("results/" + saving_folder)

            delayTargetArray = delayTargetArrayOrig * 1e-3;   % seconds

            for dIter = 1:length(delayTargetArray)
                d = delayTargetArray(dIter);

                %% Feasible PRB range from MCS limits
                [~, ~, etaMax] = parseMcs(28, 1);
                Nrb_min = ceil(pktLen / (134.16 * etaMax));

                [~, ~, etaMin] = parseMcs(0, 1);
                Nrb_max = ceil(pktLen / (134.16 * etaMin));

                NrbList = Nrb_min:1:Nrb_max;
                dvp_Harq_sim = zeros(1, length(NrbList));

                for nrbIter = 1:length(NrbList)
                    Nrb = NrbList(nrbIter);

                    %% Select minimum feasible MCS
                    etaReq = pktLen / (134 * Nrb);
                    mcs = [];
                    for mcsItem = 1:28
                        [~, ~, eta_mcsItem] = parseMcs(mcsItem, 1);
                        if eta_mcsItem >= etaReq
                            mcs = mcsItem;
                            break;
                        end
                    end

                    if isempty(mcs)
                        error("No feasible MCS found for Nrb = %d", Nrb);
                    end

                    %% PER vector from Lena-based table
                    p_Harq_all = getLenaHarqBler(simLen, mcs, SNRdb_avgRcvd, ...
                                                 hvar, M_Harq, 2, 816, Nrb);
                    p_Harq = p_Harq_all(1:M_Harq);

                    %% Queue simulation with CSO
                    [queueSize_sim, TxCount_sim, totalLatency_sim, ...
                        waitLatency_sim, numDroppedPkts] = ...
                        simulate_queue_IR_cso( ...
                            isReTxPrioritised, ...
                            f, ...
                            p_Harq, ...
                            feedbackDelay, ...
                            k_xi, ...
                            simLen, ...
                            maxQueueSize, ...
                            M_Harq, ...
                            isDetArr);

                    %% Simulation DVP
                    % Empirical ratio from returned packet latencies plus
                    % admission drops. Discarded packets have Inf latency;
                    % unresolved packets are omitted by the simulator.
                    delayTargetSlots = d / T;
                    numDelayViolations = ...
                        sum(totalLatency_sim > delayTargetSlots);
                    dvp_Harq_sim(nrbIter) = ...
                        (numDelayViolations + numDroppedPkts) / ...
                        (length(totalLatency_sim) + numDroppedPkts);
                end

                %% Save simulation results only
                writematrix(NrbList, ...
                    "results/" + saving_folder + "/" + ...
                    string(delayTargetArrayOrig(dIter)) + "NRBListSimulation.csv");

                writematrix(dvp_Harq_sim, ...
                    "results/" + saving_folder + "/" + ...
                    string(delayTargetArrayOrig(dIter)) + "DVPListSimulation.csv");
            end
        end
    end
end
