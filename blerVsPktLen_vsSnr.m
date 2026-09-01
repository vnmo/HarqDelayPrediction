clearvars -except h1 h2 h3;clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

simLen_MonteCarlo = 1e5;

Qmax = Inf;
M = Inf;
Qmax_Harq = 16;
M_Harq = 4;

T = 1e-3;
hvar = 1; % $\mu_{h^2}$
V=1;

f=1/3;

Nrb=10;
Nre = Nrb *180;
B = 180e3 * Nrb;

zeta = 1; % Decoding delay. A parameter designed later in the research.

d_Orig  = 8.5 * 1e-3;
d  = d_Orig - zeta*1e-3;
if(d<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end
delta_Orig = 2;
delta = delta_Orig + zeta;

SNRdb_avgRcvdArray = 5:5:20;
snr_avgRcvdArray = 10.^(SNRdb_avgRcvdArray/10);

pktLenArray = 8 * (20:20:250); % Packet length in bits

dvp_Arq  = zeros(1,length(pktLenArray));
dvp_Harq  = zeros(1,length(pktLenArray));

for snrIter = 1:length(snr_avgRcvdArray)
    snr_avgRcvd = snr_avgRcvdArray(snrIter)

    for pktLenIter = 1:length(pktLenArray)
        pktLen = pktLenArray(pktLenIter);

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

        kd = floor(((d/T) + delta)/(1 + delta));

        q = pArq/(1 - f);

        term1 = pArq^kd;
        term2 = f * (q)^(floor(d/T) + delta);
        term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
        term4 = 1 - f - (q)^delta;

        dvp_Arq(pktLenIter) = term1 - (term2 * term3) / term4;

        p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);
        dvp_Harq(pktLenIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);

    end
    dvp_Arq(dvp_Arq>1)=1;
    dvp_Harq(dvp_Harq>1)=1;

    figure(1)
    plot(pktLenArray/8,log10(dvp_Arq),lineOrder(3),...
        'Marker',markerOrder(1),'MarkerSize',10,...
        'color',corder(snrIter,:),'LineWidth', 3);
    hold on;
    plot(pktLenArray/8,log10(dvp_Harq),lineOrder(1),...
        'Marker',markerOrder(1),'MarkerSize',10,...
        'color',corder(snrIter,:),'LineWidth', 3);

end

figure(1)
set(gca,'Fontsize',28)
grid on;
legend( 'ARQ: SNR =  5 dB' , 'HARQ: SNR = 5 dB' ,...
        'ARQ: SNR =  10 dB' ,'HARQ: SNR = 10 dB' ,...
        'ARQ: SNR =  15 dB' ,'HARQ: SNR = 15 dB' ,...
        'ARQ: SNR =  20 dB' ,'HARQ: SNR = 20 dB' ,...
    'fontsize',32,'interpreter','latex','location','southeast','NumColumns',2,'Orientation','horizontal');
xlabel('Packet Length (B)','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');


figure(1)
saveas(gcf,'Figures/blerVsPktLen_vsSNR.fig')
saveas(gcf,'Figures/blerVsPktLen_vsSNR.epsc')
movefile('Figures/blerVsPktLen_vsSNR.epsc','Figures/blerVsPktLen_vsSNR.eps')
