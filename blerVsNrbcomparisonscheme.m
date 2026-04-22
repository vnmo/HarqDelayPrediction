close all;
clearvars -except h1 h2 h3; clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

delayTargetArrayOrig  = [2 3 5 6 8 10];
f_inverse_list = [3 6 8 10];
sinr_list = [10 20 30];
srs_list = [0];

% ===========================
% NEW: FBL-based MCS settings
% ===========================
epsilonTargetForMcs = 1e-1;   % target error rate epsilon for FBL-based MCS selection
mcsMin = 0;
mcsMax = 28;

for sinr_index = 1:length(sinr_list)
    for par_index = 1:length(f_inverse_list)
        for srs_index = 1:length(srs_list)
            
            saving_folder = "sinr"+string(sinr_list(sinr_index))+"_par"+string(f_inverse_list(par_index))+"_srs"+string(srs_list(srs_index));
            mkdir("results/"+saving_folder)

            isReTxPrioritised = true; %#ok<NASGU>
            isSimRequired = false;
            simLen = 1e6; %#ok<NASGU>
            simLen_MonteCarlo = 1e6;
            
            Qmax = Inf; %#ok<NASGU>
            M = Inf;    %#ok<NASGU>
            Qmax_Harq = 16;
            M_Harq = 4;
            
            T = 1e-3;   % slot duration
            hvar = 1;
            V = 1; %#ok<NASGU>
            
            SNRdb_avgRcvd = sinr_list(sinr_index);
            snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);

            f_inverse = f_inverse_list(par_index);
            f = 1/f_inverse;
            
            delta_Orig = 1; % Feedback delay
            pktLen = 8 * 108; % Packet length in bits

            if srs_list(srs_index) == 0
                CSO = Inf;
            else
                CSO = 1/srs_list(srs_index);
            end %#ok<NASGU>

            zeta = 1;
            delta = delta_Orig + zeta;
            delayTargetArray = (delayTargetArrayOrig) * 1e-3;

            if (min(delayTargetArray) < 1e-3)
                error("Invalid delayTarget,zeta combination. Correct it!")
            end

            % Keep same NRB sweep range logic as before
            [~, ~, etaMax] = parseMcs(28,1);
            Nrb_min = ceil(pktLen/(134.16*etaMax));
            [~, ~, etaMin] = parseMcs(0,1);
            Nrb_max = ceil(pktLen/(134.16*etaMin));
            NrbList = Nrb_min:1:Nrb_max;

            for dIter = 1:length(delayTargetArray)
                d = delayTargetArray(dIter);

                expNumTx_Harq = zeros(1,length(NrbList));
                dvp_Harq  = zeros(1,length(NrbList));
                dvp_Harq_sim = zeros(1,length(NrbList)); %#ok<NASGU>

                mcsListUsed = nan(1,length(NrbList));
                etaListUsed = nan(1,length(NrbList));
                rListUsed   = nan(1,length(NrbList));
                TtxListUsed = nan(1,length(NrbList));

                for nrbIter = 1:length(NrbList)
                    Nrb = NrbList(nrbIter);
                    Nre = Nrb * 134.16;   % channel uses available in one slot abstraction
                    B = 134e3 * Nrb;

                    % ===========================================
                    % NEW: select MCS using paper-style FBL formula
                    % ===========================================
                    [mcs, eta, r_fbl, Ttx_fbl] = selectMcsFromSinrFbl( ...
                        pktLen, snr_avgRcvd, Nre, B, epsilonTargetForMcs, mcsMin, mcsMax);

                    mcsListUsed(nrbIter) = mcs;
                    etaListUsed(nrbIter) = eta;
                    rListUsed(nrbIter)   = r_fbl;
                    TtxListUsed(nrbIter) = Ttx_fbl;

                    if isnan(mcs)
                        expNumTx_Harq(nrbIter) = Inf;
                        dvp_Harq(nrbIter) = 1;
                        continue;
                    end

                    dataBits = Nre * eta; %#ok<NASGU>

                    %% ========= HARQ =====================================
                    % KEEP THIS AS-IS: BLER from Lena
                    p_Harq = getLenaHarqBler( ...
                        simLen_MonteCarlo, mcs, SNRdb_avgRcvd, ...
                        hvar, M_Harq, 2, 816, Nrb);

                    %% Expected #transmissions
                    expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
                    for i = 1:M_Harq-1
                        expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) + ...
                            i * prod(p_Harq(1:i-1)) * (1-p_Harq(i));
                    end

                    %% DVP with HARQ
                    % KEEP FUNCTION CALL SHAPE
                    dvp_Harq(nrbIter) = getDvpHarqDM(d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, Nre);
                    % dvp_Harq(nrbIter) = getDvpHarqCSO(d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, CSO);
                end

                dvp_Harq(dvp_Harq>1) = 1;

                if ~isSimRequired
                    writematrix(NrbList,       "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"NRBListTheory.csv")
                    writematrix(dvp_Harq,      "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"DVPListTheory.csv")
                    writematrix(mcsListUsed,   "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"MCSListTheory.csv")
                    writematrix(etaListUsed,   "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"EtaListTheory.csv")
                    writematrix(rListUsed,     "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"RreqListTheory.csv")
                    writematrix(TtxListUsed,   "results/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"TtxListTheory.csv")

                    figure(1)
                    plot(NrbList, log10(dvp_Harq), lineOrder(1), ...
                        'Marker', markerOrder(dIter), 'MarkerSize', 10, ...
                        'color', corder(dIter,:), 'LineWidth', 3);
                    hold on;
                end
            end
        end
    end
end

figure(1)

%% ============================================================
% FBL-based MCS selection from SINR
% Based on the paper's normal approximation
%% ============================================================
function [mcsSel, etaSel, r_req, Ttx] = selectMcsFromSinrFbl( ...
    Lbits, snrLin, rAvail, B, epsilonTarget, mcsMin, mcsMax)

    % Shannon capacity term
    C = log2(1 + snrLin);

    % Channel dispersion
    Vdisp = (log2(exp(1)))^2 * (1 - 1/(1 + snrLin)^2);

    % Inverse Q
    qinv = qfuncinv(epsilonTarget);

    % Required channel uses from paper-style approximation
    r_req = Lbits / C + ((qinv^2) * Vdisp / (2 * C^2)) * ...
        (1 + sqrt(1 + (4 * Lbits * C) / (Vdisp * qinv^2)));

    % Corresponding transmission time T = r/B
    Ttx = r_req / B;

    % Achievable spectral efficiency using available blocklength rAvail
    % S ~= C - Q^{-1}(eps)*sqrt(V/r)
    etaFbl = C - qinv * sqrt(Vdisp / rAvail);
    etaFbl = max(0, etaFbl);

    % Choose largest MCS with eta <= etaFbl
    mcsSel = NaN;
    etaSel = NaN;

    for mcsTry = mcsMax:-1:mcsMin
        [~, ~, etaTmp] = parseMcs(mcsTry, 1);
        if etaTmp <= etaFbl
            mcsSel = mcsTry;
            etaSel = etaTmp;
            return;
        end
    end
end