clearvars -except h1 h2;clc

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
zeta = 1; % Decoding delay. A parameter designed later in the research.
% So it is not directly used. Instead d <- d-zeta, delta<- delta + zeta
delta = delta+zeta;
d=d-zeta*1e-3;
if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end
SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);

%%eta is spectral efficiency
[~, ~, etaMax] = parseMcs(28,1);
Nrb_min = ceil(pktLen/(180*etaMax));
[~, ~, etaMin] = parseMcs(0,1);
Nrb_max = ceil(pktLen/(180*etaMin));

NrbList = Nrb_min:1:Nrb_max;

Nrb = ceil((Nrb_min+Nrb_max)/2);
Nrb = 2;

NrbList = [1,2,5,15];
for NrbIter=1:length(NrbList)
 Nrb=   NrbList(NrbIter)

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
plot((fList * pktLen / T )/1000, log10(dvp_Harq),lineOrder(2), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(NrbIter,:),'LineWidth', 3);hold on
figure(1)
plot((fList * pktLen / T )/1000,rate/1000,lineOrder(2), 'Marker',markerOrder(1),'MarkerSize',10,'color',corder(NrbIter,:),'LineWidth', 3);hold on
end



% 
% 
% 
% 
% figure(1)
% set(gca,'Fontsize',28)
% grid on
% legend( '$N_{\mathrm{RB}}=1$' , '$N_{\mathrm{RB}}=2$',...
%         '$N_{\mathrm{RB}}=5$','$N_{\mathrm{RB}}=10$' ,...
%         'fontsize',32,'interpreter','latex','location','northwest','NumColumns',1,'Orientation','horizontal');
% xlabel('Arrival Rate (kbps)','fontsize',42,'interpreter','latex');
% ylabel('Throughput (kbps)','fontsize',42,'interpreter','latex');
% 
% figure(2)
% set(gca,'Fontsize',28)
% grid on
% legend( '$N_{\mathrm{RB}}=1$' , '$N_{\mathrm{RB}}=2$',...
%         '$N_{\mathrm{RB}}=5$','$N_{\mathrm{RB}}=10$' ,...
%         'fontsize',32,'interpreter','latex','location','northwest','NumColumns',1,'Orientation','horizontal');
% xlabel('Arrival Rate (kbps)','fontsize',42,'interpreter','latex');
% ylabel('DVP','fontsize',42,'interpreter','latex');
% 
% figure(1)
% saveas(gcf,'Figures/throughput_fordifferentNrb.fig')
% saveas(gcf,'Figures/throughput_fordifferentNrb.epsc')
% movefile('Figures/throughput_fordifferentNrb.epsc','Figures/throughput_fordifferentNrb.eps')
% 
% figure(2)
% saveas(gcf,'Figures/DVP_vsArrRate_vsNRB.fig')
% saveas(gcf,'Figures/DVP_vsArrRate_vsNRB.epsc')
% movefile('Figures/DVP_vsArrRate_vsNRB.epsc','Figures/DVP_vsArrRate_vsNRB.eps')