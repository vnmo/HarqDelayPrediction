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

SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);
d = 8.5e-3;
f=1/3;
delta = 2;
zeta = 1; % Decoding delay. A parameter designed later in the research.
% So it is not directly used. Instead d <- d-zeta, delta<- delta + zeta
delta = delta+zeta;
d=d-zeta*1e-3;
if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end

pktLenArray = 8 * [30 50 100 200]; % Packet length in bits
for pktLenIter = 1:length(pktLenArray)
    pktLen = pktLenArray(pktLenIter)

    % eta is spectral efficiency
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

    for nrbIter = 1:length(NrbList)
        Nrb=NrbList(nrbIter);
        Nre = Nrb *180;
        B = 180e3 * Nrb;

        % Select the minimum possible MCS to transfer the bits
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
        %         p(nrbIter)= 1-exp((1-exp(dataBits*log(2)/(T*B)))/(snr*var_hSquare));

        % FBL Approximation
        pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);

        %% Expected #transmissions
        expNumTx_Arq(nrbIter) = 1/(1-(pArq));


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
            [~,~,latencyInSlotsVector,~] = simulate_queue_noIR(isReTxPrioritised,f,pArq(nrbIter),delta,simLen,Qmax,M,false);
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
        plot(NrbList,log10(dvp_Arq),lineOrder(pktLenIter),'Marker',markerOrder(pktLenIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
        hold on
        plot(NrbList,log10(dvp_Arq_sim),lineOrder(pktLenIter),'Marker',markerOrder(pktLenIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);

        figure(2)
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(pktLenIter),'Marker',markerOrder(pktLenIter),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
        hold on
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq_sim),lineOrder(pktLenIter),'Marker',markerOrder(pktLenIter),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);
    else
        figure(1)
        plot(NrbList,log10(dvp_Arq),lineOrder(3),...
            'Marker',markerOrder(pktLenIter),'MarkerSize',10,...
            'color',corder(pktLenIter,:),'LineWidth', 3);
        hold on;
        plot(NrbList,log10(dvp_Harq),lineOrder(1),...
            'Marker',markerOrder(pktLenIter),'MarkerSize',10,...
            'color',corder(pktLenIter,:),'LineWidth', 3);

        figure(2)
        plot(NrbList.*expNumTx_Arq,log10(dvp_Arq),lineOrder(3),...
            'Marker',markerOrder(pktLenIter),'MarkerSize',10,...
            'color',corder(pktLenIter,:),'LineWidth', 3);
        hold on;
        plot(NrbList.*expNumTx_Harq,log10(dvp_Harq),lineOrder(1),...
            'Marker',markerOrder(pktLenIter),'MarkerSize',10,...
            'color',corder(pktLenIter,:),'LineWidth', 3);

    end

end

figure(1)
set(gca,'Fontsize',28)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% title(titletext,'fontsize',42,'interpreter','latex');
grid on;
legend( 'ARQ: $n = 30$ B' ,'HARQ: $n = 30$ B' ,...
    'ARQ: $n = 50$ B' ,'HARQ: $n = 50$ B' ,...
    'ARQ: $n = 100$ B','HARQ: $n = 100$ B',...
    'ARQ: $n = 200$ B','HARQ: $n = 200$ B',...
    'fontsize',32,'interpreter','latex','location','northeast','NumColumns',2,'Orientation','horizontal');

xlabel('Allocated $N_{RB}$ per slot','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
% grid minor
ylim([-4,0])
xlim([0,40])

figure(2)
set(gca,'Fontsize',28)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. average consumed resources per packet.";
% title(titletext,'fontsize',42,'interpreter','latex');
grid on;
legend( 'ARQ: $n = 30$ B' ,'HARQ: $n = 30$ B' ,...
    'ARQ: $n = 50$ B' ,'HARQ: $n = 50$ B' ,...
    'ARQ: $n = 100$ B','HARQ: $n = 100$ B',...
    'ARQ: $n = 200$ B','HARQ: $n = 200$ B',...
    'fontsize',32,'interpreter','latex','location','northeast','NumColumns',2,'Orientation','horizontal');

xlabel('Average RB consumed per packet','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
% grid minor
ylim([-4,0])
xlim([0,40])


figure(1)
saveas(gcf,'Figures/blerVsNrb_vsPktLen.fig')
saveas(gcf,'Figures/blerVsNrb_vsPktLen.epsc')
movefile('Figures/blerVsNrb_vsPktLen.epsc','Figures/blerVsNrb_vsPktLen.eps')

figure(2)
saveas(gcf,'Figures/blerVsNrbUsed_vsPktLen.fig')
saveas(gcf,'Figures/blerVsNrbUsed_vsPktLen.epsc')
movefile('Figures/blerVsNrbUsed_vsPktLen.epsc','Figures/blerVsNrbUsed_vsPktLen.eps')