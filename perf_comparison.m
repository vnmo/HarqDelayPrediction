clearvars -except h1;clc

isReTxPrioritised = true;

simLen =1e6;
simLen_MonteCarlo = 1e5;

Qmax = Inf;
M = Inf;
Qmax_Harq = 16;
M_Harq = 4;

T = 1e-3;
hvar = 1;
V=1;

pktLen = 8 * 100; % Packet length in bits
d = 8.5e-3;
f=1/3;


delta = 2;
zeta = 1; % Decoding delay. A parameter designed later in the research.
delta = delta+zeta;
d_orig=d;
d = d - zeta * T;
if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end


SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);

[~, ~, etaMax] = parseMcs(28,1);
Nrb_min = ceil(pktLen/(180*etaMax));
[~, ~, etaMin] = parseMcs(0,1);
Nrb_max = ceil(pktLen/(180*etaMin));

NrbList = Nrb_min:1:Nrb_max;

dvp_Arq  = zeros(1,length(NrbList));
dvp_Arq_sim = zeros(1,length(NrbList));
dvp_Harq  = zeros(1,length(NrbList));
dvp_Harq_sim = zeros(1,length(NrbList));

dvp_Arq_IF =zeros(1,length(NrbList));
dvp_Harq_IF =zeros(1,length(NrbList));
dvp_Arq_TAF =zeros(1,length(NrbList));


for nrbIter = 1:length(NrbList)
    Nrb=NrbList(nrbIter);
    Nre = Nrb *180;
    B = 180e3 * Nrb;

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

    pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);
    p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);


    kd = floor(((d/T) + delta)/(1 + delta));
    q = pArq/(1 - f); % Let

    term1 = pArq^kd;
    term2 = f * (q)^(floor(d/T) + delta);
    term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
    term4 = 1 - f - (q)^delta;

    dvp_Arq(nrbIter) = term1 - (term2 * term3) / term4;

    if (f > 1-pArq)
        dvp_Arq_sim(nrbIter) = 1;
    else
        [~,~,latencyInSlotsVector,~] = simulate_queue_noIR(isReTxPrioritised,f,pArq,delta,simLen,Qmax,M,false);
        delaytargetInSlots = (d/T);
        dvp_Arq_sim(nrbIter) = (sum(latencyInSlotsVector > delaytargetInSlots))/length(latencyInSlotsVector);
    end

    delta_temp = 0; % Make RTT = 1;
    
    kd = floor(((d_orig/T) + delta_temp)/(1 + delta_temp));
    q = pArq/(1 - f); % Let

    term1 = pArq^kd;
    term2 = f * (q)^(floor(d_orig/T) + delta_temp);
    term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta_temp));
    term4 = 1 - f - (q)^delta_temp;

    dvp_Arq_IF(nrbIter) = term1 - (term2 * term3) / term4;


    dvp_Harq(nrbIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);

    [~,~,latencyInSlotsVector,~] = simulate_queue_IR(isReTxPrioritised,f,p_Harq,delta,simLen,Qmax_Harq,M_Harq,false);

    delaytargetInSlots = (d/T);
    dvp_Harq_sim(nrbIter) = (sum(latencyInSlotsVector > delaytargetInSlots))/length(latencyInSlotsVector);

    delta_temp = 0; % Make RTT = 1;
    dvp_Harq_IF(nrbIter) = getDvpHarq (d_orig, delta_temp, p_Harq, f, T, M_Harq, Qmax_Harq);








end

dvp_Arq(dvp_Arq>1)=1;
dvp_Arq_sim(dvp_Arq_sim>1)=1;
dvp_Harq(dvp_Harq>1)=1;
dvp_Harq_sim(dvp_Harq_sim>1)=1;
dvp_Arq_IF(dvp_Arq_IF>1)=1;
dvp_Harq_IF(dvp_Harq_IF>1)=1;



corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

plot(NrbList,log10(dvp_Harq),    lineOrder(1), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);hold on;
plot(NrbList,log10(dvp_Arq),     lineOrder(3), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3)
plot(NrbList,log10(dvp_Harq_IF), lineOrder(1), 'Marker',markerOrder(6),'MarkerSize',10,'color',corder(6,:),'LineWidth', 3)
plot(NrbList,log10(dvp_Arq_IF),  lineOrder(3), 'Marker',markerOrder(4),'MarkerSize',10,'color',corder(4,:),'LineWidth', 3)
plot(NrbList,log10(dvp_Harq_sim),lineOrder(2), 'Marker',markerOrder(2),'MarkerSize',10,'color',corder(5,:),'LineWidth', 3)
plot(NrbList,log10(dvp_Arq_sim), lineOrder(2), 'Marker',markerOrder(3),'MarkerSize',10,'color',corder(3,:),'LineWidth', 3)

set(gca,'Fontsize',28)
grid on
legend( 'HARQ' , 'ARQ',...
        'HARQ-IF','ARQ-IF' ,...
        'HARQ: Simulation','ARQ: Simulation',...
        'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');
xlabel('Allocated $N_{RB}$ per slot','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
ylim([-12,0])
xlim([0,20])

figure(1)
saveas(gcf,'Figures/comparison.fig')
saveas(gcf,'Figures/comparison.epsc')
movefile('Figures/comparison.epsc','Figures/comparison.eps')

