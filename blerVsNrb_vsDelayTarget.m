close all;
clearvars -except h1 h2 h3;clc

lineOrder = ["-","--",":","-","--",":","-.","-","--",":","-","--",":","-.",];
markerOrder = ["none","o","x","diamond","|","square","none","o","x","diamond","|","square",];

delayTargetArrayOrig  = [2 3 5 6 7 8 10 11 12];
f_inverse_list = [3 6 8 10];
sinr_list = [10 20];
srs_list = [80 320];

for sinr_index = 1:length(sinr_list)
    for par_index = 1:length(f_inverse_list)
        for srs_index = 1:length(srs_list)
            
            saving_folder = "sinr"+string(sinr_list(sinr_index))+"_par"+string(f_inverse_list(par_index))+"_srs"+string(srs_list(srs_index));
            mkdir("results_phase_fifo/"+saving_folder)
            isReTxPrioritised = true;
            isSimRequired = false;
            simLen =1e7;
            simLen_MonteCarlo = 1e7;
            
            Qmax = Inf;
            M = Inf;
            Qmax_Harq = 16;
            M_Harq = 4;
            
            T = 1e-3;
            hvar = 1; % $\mu_{h^2}$
            V=1;
            
            SNRdb_avgRcvd = sinr_list(sinr_index);

            snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);
            f_inverse = f_inverse_list(par_index);
            f=1/f_inverse;
            
            delta_Orig = 1.0; % Feedback delay
            pktLen = 8 * 108; % Packet length in bits
            CSO = 1/srs_list(srs_index); % Control Signalling Overhead
            
            
            zeta = 1.0; % Raw decoding delay in slots, rounded up in service analysis
            delta = 1; % Raw feedback delay in slots, rounded up separately
            delayTargetArray=(delayTargetArrayOrig)*1e-3;
            if(min(delayTargetArray)<1e-3)
                error("Invalid delayTarget,zeta combination. Correct it!")
            end
            for dIter = 1:length(delayTargetArray)
                d = delayTargetArray(dIter)
            
                [~, ~, etaMax] = parseMcs(28,1);
                Nrb_min = ceil(pktLen/(134.16*etaMax));
                [~, ~, etaMin] = parseMcs(0,1);
                Nrb_max = ceil(pktLen/(134.16*etaMin));
            
                NrbList = Nrb_min:1:Nrb_max;
            
                expNumTx_Arq = zeros(1,length(NrbList));
                dvp_Arq  = zeros(1,length(NrbList));
                dvp_Arq_sim = zeros(1,length(NrbList));
                expNumTx_Harq = zeros(1,length(NrbList));
                dvp_Harq  = zeros(1,length(NrbList));
                dvp_Harq_sim = zeros(1,length(NrbList));
                McsPrev = -1;
                for nrbIter = 1:length(NrbList)
                    Nrb=NrbList(nrbIter);
                    Nre = Nrb * 134.16;
                    B = 134e3 * Nrb;
            
                    etaMin = pktLen/(134*Nrb);
                    eta=[];
                    mcs=[];
                    for mcsItem =1:28
                        [~, ~, eta_mcsItem] = parseMcs(mcsItem,1);
                        if (eta_mcsItem>=etaMin)
                            eta = eta_mcsItem;
                            mcs=mcsItem;
                            break;
                        end
                    end
                    mcs
                    McsPrev = mcs;
                    dataBits = Nre *eta;
            
                    
            
                    p_Harq = getLenaHarqBler(simLen_MonteCarlo, mcs, SNRdb_avgRcvd, hvar, M_Harq, 2, 816, Nrb);
                    expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
                    for i = 1:M_Harq-1
                        expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) + i * prod(p_Harq(1:i-1)) * (1-p_Harq(i));
                    end
                    dvp_Harq(nrbIter) = getDvpHarqCSO(d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, CSO, zeta);
            
                end
                dvp_Harq(dvp_Harq>1)=1;
            
                if(isSimRequired)
            
                else
                    
                    writematrix(NrbList, "results_phase_fifo/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"NRBListTheory.csv")
                    writematrix(dvp_Harq, "results_phase_fifo/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"DVPListTheory.csv")
                    

                    
            
                    plot(NrbList,log10(dvp_Harq),lineOrder(1),...
                        'Marker',markerOrder(dIter),'MarkerSize',10,...
                        'LineWidth', 3);
                    hold on;
            
                end
            
            end
        end
    end
end


figure(1)
