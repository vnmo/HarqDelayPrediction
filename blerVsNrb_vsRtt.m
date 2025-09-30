% close all;
clearvars -except h1 h2 h3;clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

isReTxPrioritised = true;
isSimRequired = false;
simLen =1e4;
simLen_MonteCarlo = 1e5;

Qmax = Inf;
M = Inf;
Qmax_Harq = 16;
M_Harq = 4;

T = 1e-3;
hvar = 1; % $\mu_{h^2}$
V=1;

pktLen = 8 * 100; % Packet length in bits
SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);
d = 8.5e-3;
f=1/3;

%%eta is spectral efficiency
[~, ~, etaMax] = parseMcs(28,1);
Nrb_min = ceil(pktLen/(180*etaMax));
[~, ~, etaMin] = parseMcs(0,1);
Nrb_max = ceil(pktLen/(180*etaMin));
NrbList = Nrb_min:1:Nrb_max;

expNumTx_Arq = zeros(1,length(NrbList));
dvp_Arq  = zeros(1,length(NrbList));
dvp_Arq_sim = zeros(1,length(NrbList));
expNumTx_Harq = zeros(1,length(NrbList));
dvp_Harq  = zeros(1,length(NrbList));
dvp_Harq_sim = zeros(1,length(NrbList));


deltaArray = [0 2 4];

zeta = 0; % Decoding delay. A parameter designed later in the research.
% So it is not directly used. Instead d <- d-zeta, delta<- delta + zeta
deltaArray = deltaArray + zeta;
rttArray = 1 + deltaArray;
d = d - zeta * T;

if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end
for rttIter = 1:length(rttArray)
    delta = deltaArray(rttIter);
    rtt = rttArray(rttIter)

    for nrbIter = 1:length(NrbList)
        Nrb=NrbList(nrbIter);
        Nre = Nrb *180;
        B = 180e3 * Nrb;

        %Select the minimum possible MCS to transfer the bits
        etaMin = pktLen/(180*Nrb);
        eta=[];
        mcs=[];
        for mcsItem =0:28
            [~, ~, eta_mcsItem] = parseMcs(mcsItem,1);
            if (eta_mcsItem>=etaMin)
                eta = eta_mcsItem;
                mcs=mcsItem;
                break;
            end
        end

        dataBits = Nre *eta;

        %% =========ARQ========================================================

        %% PER computation.
        %         % Shannon Outage
        %         p= 1-exp((1-exp(dataBits*log(2)/(T*B)))/(snr*var_hSquare));

        % FBL Approximation
        pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);

        %% Expected #transmissions
        expNumTx_Arq(nrbIter) = 1/(1-pArq);


        %% ARQ DVP with feedback delay (Theory)
        % total delay closed form
        kd = floor(((d/T) + delta)/(1 + delta));

        % Let
        q = pArq/(1 - f);

        term1 = pArq^kd;
        term2 = f * (q)^(floor(d/T) + delta);
        term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
        term4 = 1 - f - (q)^delta;

        dvp_Arq(nrbIter) = term1 - (term2 * term3) / term4;


        %% ARQ DVP with feedback delay (Simulation)
        if(isSimRequired)
            [~,~,latencyInSlotsVector,~] = simulate_queue_noIR(isReTxPrioritised,f,pArq,delta,simLen,Qmax,M,false);
            delaytargetInSlots = (d/T);
            dvp_Arq_sim(nrbIter) = (sum(latencyInSlotsVector > delaytargetInSlots))/length(latencyInSlotsVector);
        end

        %% =========HARQ========================================================

        %% PER computation.
        p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);

        %% Expected #transmissions
        expNumTx_Harq(nrbIter) = M_Harq * prod(p_Harq);
        for i = 1:M_Harq-1
            expNumTx_Harq(nrbIter) = expNumTx_Harq(nrbIter) + i * prod(p_Harq(1:i-1)) * (1-p_Harq(i));
        end

        %% DVP of with HARQ; (Truncated at M = 4)
        dvp_Harq(nrbIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);
    end

    dvp_Arq(dvp_Arq>1)=1;
    dvp_Harq(dvp_Harq>1)=1;

    if(isSimRequired)
        figure(1)
        plot(NrbList,log10(dvp_Arq),lineOrder(rttIter),'Marker',markerOrder(rttIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
        hold on
        plot(NrbList,log10(dvp_Arq_sim),lineOrder(rttIter),'Marker',markerOrder(rttIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);

        figure(2)
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(rttIter),'Marker',markerOrder(rttIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
        hold on
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq_sim),lineOrder(rttIter),'Marker',markerOrder(rttIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);
    else
        figure(1)
        plot(NrbList,log10(dvp_Arq),lineOrder(3),...
            'Marker',markerOrder(rttIter),'MarkerSize',10,...
            'color',corder(rttIter,:),'LineWidth', 3);
        hold on;
        plot(NrbList,log10(dvp_Harq),lineOrder(1),...
            'Marker',markerOrder(rttIter),'MarkerSize',10,...
            'color',corder(rttIter,:),'LineWidth', 3);

        figure(2)
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(3),...
            'Marker',markerOrder(rttIter),'MarkerSize',10,...
            'color',corder(rttIter,:),'LineWidth', 3);
        hold on;
        plot(NrbList.*expNumTx_Harq,log10(dvp_Harq),lineOrder(1),...
            'Marker',markerOrder(rttIter),'MarkerSize',10,...
            'color',corder(rttIter,:),'LineWidth', 3);

    end

end

figure(1)
set(gca,'Fontsize',28)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% title(titletext,'fontsize',42,'interpreter','latex');
grid on;
legend( 'ARQ: RTT = 1','HARQ: RTT =  1',...
        'ARQ: RTT = 3','HARQ: RTT =  3',...
        'ARQ: RTT = 5','HARQ: RTT =  5',...
    'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');

xlabel('Allocated $N_{RB}$ per slot','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
% grid minor
ylim([-12,0])
xlim([0,20])


figure(2)
set(gca,'Fontsize',28)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. average consumed resources per packet.";
% title(titletext,'fontsize',42,'interpreter','latex');
grid on;
legend( 'ARQ: RTT = 1 (IF)','HARQ: RTT =  1 (IF)',...
        'ARQ: RTT = 3','HARQ: RTT =  3',...
        'ARQ: RTT = 5','HARQ: RTT =  5',...
    'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');

xlabel('Average $N_{RB}$ per packet','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
% grid minor
ylim([-12,0])
xlim([0,20])


figure(1)
saveas(gcf,'Figures/blerVsNrb_vsRtt.fig')
saveas(gcf,'Figures/blerVsNrb_vsRtt.epsc')
movefile('Figures/blerVsNrb_vsRtt.epsc','Figures/blerVsNrb_vsRtt.eps')


figure(2)
saveas(gcf,'Figures/blerVsNrbUsed_vsRtt.fig')
saveas(gcf,'Figures/blerVsNrbUsed_vsRtt.epsc')
movefile('Figures/blerVsNrbUsed_vsRtt.epsc','Figures/blerVsNrbUsed_vsRtt.eps')