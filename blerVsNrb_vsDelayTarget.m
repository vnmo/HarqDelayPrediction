% Proposed model: analytical DVP sweeps using getDvpHarqCSO.
% Outputs are *Theory.csv; this script does not run a queueing simulation.
% The results_phase_fifo root keeps revised curves apart from old outputs.
% The corresponding simulation driver is blerVsNrb_vsDelayTarget_queuingsim.
% getLenaHarqBler supplies link-level PER inputs, not queueing observations.
% Scenario settings are defined below; see CSO_REVISION_NOTES.md for roles.
close all;
clearvars -except h1 h2 h3;clc

% corder = [colororder colororder];
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
            % This flag controls the existing output branch only; setting
            % it true does not call simulate_queue_IR_cso in this script.
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
            % delayTargetArrayOrig  = [2.5 5.5 8.5 11.5 14.5] * 1e-3;
            
            
            zeta = 1.0; % Raw decoding delay in slots, rounded up in service analysis
            delta = 1; % Raw feedback delay in slots, rounded up separately
            delayTargetArray=(delayTargetArrayOrig)*1e-3;
            if(min(delayTargetArray)<1e-3)
                error("Invalid delayTarget,zeta combination. Correct it!")
            end
            for dIter = 1:length(delayTargetArray)
                d = delayTargetArray(dIter)
            
                % eta is spectral efficiency
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
            
                    % Select the minimum possible MCS to transfer the bits
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
            
                    
                    %% =========HARQ========================================================
            
                    %% PER computation.
                    % p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);
                    p_Harq = getLenaHarqBler(simLen_MonteCarlo, mcs, SNRdb_avgRcvd, hvar, M_Harq, 2, 816, Nrb);
                    % p_Harq = [0.999, 0.999, 0.999, 0.99];
                    %% Expected #transmissions
                    expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
                    for i = 1:M_Harq-1
                        expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) + i * prod(p_Harq(1:i-1)) * (1-p_Harq(i));
                    end
                    %% Proposed DVP approximation with the configured attempt cap
                    % dvp_Harq(nrbIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);
                    dvp_Harq(nrbIter) = getDvpHarqCSO(d, delta, p_Harq, f, T, M_Harq, Qmax_Harq, CSO, zeta);
            
                end
                % dvp_Harq
                % dvp_Arq(dvp_Arq>1)=1;
                dvp_Harq(dvp_Harq>1)=1;
            
                if(isSimRequired)
                    % figure(1)
                    % plot(NrbList,log10(dvp_Arq),lineOrder(dIter),'Marker',markerOrder(dIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
                    % hold on
                    % plot(NrbList,log10(dvp_Arq_sim),lineOrder(dIter),'Marker',markerOrder(dIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);
            
                    % figure(2)
                    % plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(dIter),'Marker',markerOrder(dIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
                    % hold on
                    % plot(NrbList.*expNumTx_Arq,log10(dvp_Arq_sim),lineOrder(dIter),'Marker',markerOrder(dIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);
                else
                    % figure(1)
                    % plot(NrbList,log10(dvp_Arq),lineOrder(3),...
                    %     'Marker',markerOrder(dIter),'MarkerSize',10,...
                    %     'color',corder(dIter,:),'LineWidth', 3);
                    % hold on;
                    
                    writematrix(NrbList, "results_phase_fifo/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"NRBListTheory.csv")
                    writematrix(dvp_Harq, "results_phase_fifo/"+saving_folder+"/"+string(delayTargetArrayOrig(dIter))+"DVPListTheory.csv")
                    

                    
            
                    plot(NrbList,log10(dvp_Harq),lineOrder(1),...
                        'Marker',markerOrder(dIter),'MarkerSize',10,...
                        'LineWidth', 3);
                    hold on;
                    % figure(2)
                    % plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(3),...
                    %     'Marker',markerOrder(dIter),'MarkerSize',10,...
                    %     'color',corder(dIter,:),'LineWidth', 3);
                    % hold on;
                    % plot(NrbList.*expNumTx_Harq,log10(dvp_Harq),lineOrder(1),...
                    %     'Marker',markerOrder(dIter),'MarkerSize',10,...
                    %     'color',corder(dIter,:),'LineWidth', 3);
            
                end
            
            end
        end
    end
end
% figure(1)
% set(gca,'Fontsize',28)
% % titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% % title(titletext,'fontsize',36,'interpreter','latex');
% grid on;
% legend( 'HARQ: $d = 2$ ms',...
%     'HARQ: $d = 3$ ms' ,...
%     'HARQ: $d = 5$ ms' ,...
%     'HARQ: $d = 6$ ms',...
%     'HARQ: $d = 14.5$ ms',...
%     'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');
% 
% xlabel('PRB per slot','fontsize',42,'interpreter','latex');
% ylabel('DVP','fontsize',42,'interpreter','latex');
% % grid minor
% ylim([-6,0])
% xlim([0,20])
% yticks([-6,-5,-4,-3,-2,-1,0])
% yticklabels({'10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}','1'})


figure(1)
% saveas(gcf,'Figures/blerVsNrb_vsDelayTarget.fig')
% saveas(gcf,'Figures/blerVsNrb_vsDelayTarget.epsc')
% movefile('Figures/blerVsNrb_vsDelayTarget.epsc','Figures/blerVsNrb_vsDelayTarget.eps')
% 
% figure(2)
% saveas(gcf,'Figures/blerVsNrbUsed_vsDelayTarget.fig')
% saveas(gcf,'Figures/blerVsNrbUsed_vsDelayTarget.epsc')
% movefile('Figures/blerVsNrbUsed_vsDelayTarget.epsc','Figures/blerVsNrbUsed_vsDelayTarget.eps')
