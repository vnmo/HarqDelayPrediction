hold off;clear;clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

simLen_MonteCarlo = 1e5;

isReTxPrioritised = true;
maxQueueSize = Inf;
MaxTxPerPkt = Inf;

T = 1e-3;
hvar = 10; % $\mu_{h^2}$
V=1;

pkTLen = 8 * 1000; % Packet length in bits

%%eta is spectral efficiency
[~, ~, etaMax] = parseMcs(28,1);
Nrb_min = ceil(pkTLen/(180*etaMax));
[~, ~, etaMin] = parseMcs(0,1);
Nrb_max = ceil(pkTLen/(180*etaMin));

NrbList = Nrb_min:1:Nrb_max;

save('runTimeData.mat')

RcvdSnrAvgDbList = [10 20];
for SNRdbIter = 1:length(RcvdSnrAvgDbList)
    SNRdb = RcvdSnrAvgDbList(SNRdbIter);
    snrRcvd_avg = 10.^(SNRdb/10);

for nrbIter = 1:length(NrbList)
    Nrb=NrbList(nrbIter);
    Nre = Nrb *180;
    B = 180e3 * Nrb;

    %Select the minimum possible MCS to transfer the bits
    etaMin = pkTLen/(180*Nrb);
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



    %% PER computation.
    snrTx_avg = snrRcvd_avg/hvar;
    % Shannon Outage
    p_shannon(nrbIter)= 1-exp((1-exp(dataBits*log(2)/(T*B)))/(snrTx_avg*hvar));

    % FBL Approximation
    p_FBL(nrbIter) = getFblPer(simLen_MonteCarlo,snrRcvd_avg,eta,Nre,V,hvar);

end


figure(1)
plot(NrbList,log10(p_shannon),lineOrder(1),...
    'Marker',markerOrder(2),'MarkerSize',12,...
    'color',corder(1,:),'LineWidth', 3);
hold on
plot(NrbList,log10(p_FBL),lineOrder(1),...
    'Marker',markerOrder(1),'MarkerSize',12,...
    'color',corder(2,:),'LineWidth', 3);
end


set(gca,'Fontsize',24)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% title(titletext,'fontsize',36,'interpreter','latex');
grid on;
legend('SNR = 10dB, Shannon outage ','SNR = 10dB, FBL approximation',...
    'SNR = 20dB, Shannon outage ','SNR = 20dB, FBL approximation',...
    'fontsize',24,'interpreter','latex','location','northeast','NumColumns',2,'Orientation','horizontal');

xlabel('$N_{RB}$','fontsize',36,'interpreter','latex');
ylabel('log$_{10}(p)$','fontsize',36,'interpreter','latex');
grid minor















clear;clc
load('runTimeData.mat')
RcvdSnrAvgDbList = 0:3:30;
NrbList = NrbList(floor([1 2]*length(NrbList)/3));

for nRbIter = 1:length(NrbList)
    Nrb=NrbList(nRbIter);
    Nre = Nrb *180;
    B = 180e3 * Nrb;

    % Select the minimum possible MCS to transfer the bits
    etaMin = pkTLen/(180*Nrb);
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

for snrIter = 1:length(RcvdSnrAvgDbList)
    SNRdb = RcvdSnrAvgDbList(snrIter);
    snrRcvd_avg = 10.^(SNRdb/10);

    %% PER computation.
    snrTx_avg = snrRcvd_avg/hvar;
    % Shannon Outage
    p_shannon(snrIter)= 1-exp((1-exp(dataBits*log(2)/(T*B)))/(snrTx_avg*hvar));

    % FBL Approximation
    p_FBL(snrIter) = getFblPer(simLen_MonteCarlo,snrRcvd_avg,eta,Nre,V,hvar);

end


figure(2)
plot(RcvdSnrAvgDbList,log10(p_shannon),lineOrder(1),...
    'Marker',markerOrder(2),'MarkerSize',12,...
    'color',corder(1,:),'LineWidth', 3);
hold on
plot(RcvdSnrAvgDbList,log10(p_FBL),lineOrder(2),...
    'Marker',markerOrder(1),'MarkerSize',12,...
    'color',corder(2,:),'LineWidth', 3);
end


set(gca,'Fontsize',24)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% title(titletext,'fontsize',36,'interpreter','latex');
grid on;
legend("$N_{\mathrm{RB}}="+ NrbList(1)+"$, Shannon outage ","$N_{\mathrm{RB}}="+ NrbList(1)+"$, FBL approximation",...
    "$N_{\mathrm{RB}}="+ NrbList(2)+"$, Shannon outage ","$N_{\mathrm{RB}}="+ NrbList(2)+"$, FBL approximation",...
    'fontsize',24,'interpreter','latex','location','northeast','NumColumns',2,'Orientation','horizontal');

xlabel('SNR (dB)','fontsize',36,'interpreter','latex');
ylabel('log$_{10}(p)$','fontsize',36,'interpreter','latex');
grid minor




figure(1)
saveas(gcf,'Figures/Validation_FBL_1.fig')
saveas(gcf,'Figures/Validation_FBL_1.epsc')

figure(2)
saveas(gcf,'Figures/Validation_FBL_2.fig')
saveas(gcf,'Figures/Validation_FBL_2.epsc')