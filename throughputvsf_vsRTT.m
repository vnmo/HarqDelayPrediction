clearvars -except h1 h2 h3;clc

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
delta = 2;
zeta = 0; % Decoding delay. A parameter designed later in the research.
% So it is not directly used. Instead d <- d-zeta, delta<- delta + zeta
deltaArray = [0 2 4 6];
deltaArray = deltaArray + zeta;
d=d-zeta*1e-3;
if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end
rttArray = 1+deltaArray;


SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);

%%eta is spectral efficiency
[~, ~, etaMax] = parseMcs(28,1);
Nrb_min = ceil(pktLen/(180*etaMax));
[~, ~, etaMin] = parseMcs(0,1);
Nrb_max = ceil(pktLen/(180*etaMin));

NrbList = Nrb_min:1:Nrb_max;


Nrb=   10;
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


for rttIter=1:length(rttArray)

    rtt =rttArray(rttIter)
    delta= deltaArray(rttIter);


    %% ========= PER COMPUTATION =========
    pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);
    p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);


    fList= 0.01:0.01:0.99;

    dvp_Arq  = zeros(1,length(fList));
    dvp_Harq  = zeros(1,length(fList));

    for fIter=1:length(fList)
        f=fList(fIter);
        % %% ========= ARQ =========
        % kd = floor(((d/T) + delta)/(1 + delta));
        % q = pArq/(1 - f); % Let
        %
        % term1 = pArq^kd;
        % term2 = f * (q)^(floor(d/T) + delta);
        % term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
        % term4 = 1 - f - (q)^delta;
        %
        % dvp_Arq(fIter) = term1 - (term2 * term3) / term4;


        %% ========= HARQ =========
        dvp_Harq(fIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);

    end

    % dvp_Arq(dvp_Arq>1)=1;
    dvp_Harq(dvp_Harq>1)=1;
    rate = (fList * pktLen / T ) .* (1-dvp_Harq);
   
    corder = colororder;
    lineOrder = ["-","--",":","-","--",":","-."];
    markerOrder = ["none","o","x","diamond","|","square"];


    figure(2)
    plot((fList * pktLen / T )/1000, log10(dvp_Harq),lineOrder(1), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(rttIter,:),'LineWidth', 3);hold on
    figure(1)
    plot((fList * pktLen / T )/1000,rate/1000,lineOrder(1), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(rttIter,:),'LineWidth', 3);hold on

end

figure(1)

set(gca,'Fontsize',28)
grid on
legend('$N_{\mathrm{RB}}=1,\;\, RTT=4$' , '$N_{\mathrm{RB}}=2,\;\, RTT=4$','$N_{\mathrm{RB}}=5,\;\, RTT=4$','$N_{\mathrm{RB}}=15, RTT=4$',...
       '$N_{\mathrm{RB}}=10, RTT=1$' , '$N_{\mathrm{RB}}=10, RTT=3$','$N_{\mathrm{RB}}=10, RTT=5$','$N_{\mathrm{RB}}=10, RTT=7$' ,...
'fontsize',32,'interpreter','latex','location','northwest','NumColumns',2,'Orientation','vertical');
xlabel('Arrival Rate (kbps)','fontsize',42,'interpreter','latex');
ylabel('Throughput (kbps)','fontsize',42,'interpreter','latex');


figure(2)
set(gca,'Fontsize',28)
grid on
legend('$N_{\mathrm{RB}}=1,\;\, RTT=4$' , '$N_{\mathrm{RB}}=2,\;\, RTT=4$','$N_{\mathrm{RB}}=5,\;\, RTT=4$','$N_{\mathrm{RB}}=15, RTT=4$',...
       '$N_{\mathrm{RB}}=10, RTT=1$' , '$N_{\mathrm{RB}}=10, RTT=3$','$N_{\mathrm{RB}}=10, RTT=5$','$N_{\mathrm{RB}}=10, RTT=7$' ,...
'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','vertical');
xlabel('Arrival Rate (kbps)','fontsize',42,'interpreter','latex');
ylabel('Throughput (kbps)','fontsize',42,'interpreter','latex');

figure(1)
saveas(gcf,'Figures/throughput_fordifferentRTTAndNrb.fig')
saveas(gcf,'Figures/throughput_fordifferentRTTAndNrb.epsc')
movefile('Figures/throughput_fordifferentRTTAndNrb.epsc','Figures/throughput_fordifferentRTTAndNrb.eps')

figure(2)
saveas(gcf,'Figures/DVP_vsArrRate_vsRTTAndNrb.fig')
saveas(gcf,'Figures/DVP_vsArrRate_vsRTTAndNrb.epsc')
movefile('Figures/DVP_vsArrRate_vsRTTAndNrb.epsc','Figures/DVP_vsArrRate_vsRTTAndNrb.eps')