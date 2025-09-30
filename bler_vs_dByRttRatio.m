% close all;
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

SNRdb_avgRcvd = 10;
snr_avgRcvd = 10.^(SNRdb_avgRcvd/10);
f=1/3;
pktLen = 8 * 100; % Packet length in bits

Nrb=10;
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

zeta = 0; % Decoding delay. A parameter designed later in the research.
% So it is not directly used. Instead d <- d-zeta, delta<- delta + zeta

delayTargetArray  = (2.5:0.25:20) * 1e-3;
delayTargetArray  = delayTargetArray - zeta*1e-3;
if(min(delayTargetArray)<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end
deltaArray = 0:2:6;
deltaArray = deltaArray + zeta;
rttArray = 1 + deltaArray;

for rttIter = 1:length(rttArray)
    delta = deltaArray(rttIter);
    rtt = rttArray(rttIter)

    dvp_Arq  = zeros(1,length(delayTargetArray));
    dvp_Harq  = zeros(1,length(delayTargetArray));


    for dIter = 1:length(delayTargetArray)
        d = delayTargetArray(dIter);

        %% =========ARQ========================================================
        % PER computation.
        %         % Shannon Outage
        %         p(nrbIter)= 1-exp((1-exp(dataBits*log(2)/(T*B)))/(snr*var_hSquare));
        % FBL Approximation
        pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);

        % ARQ DVP with feedback delay (Theory)
        % total delay closed form
        kd = floor(((d/T) + delta)/(1 + delta));

        % Let
        q = pArq/(1 - f);

        term1 = pArq^kd;
        term2 = f * (q)^(floor(d/T) + delta);
        term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
        term4 = 1 - f - (q)^delta;

        dvp_Arq(dIter) = term1 - (term2 * term3) / term4;

        %% =========HARQ========================================================
        % PER computation.
        p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);
        % DVP of with HARQ; (Truncated at M = 4)
        dvp_Harq(dIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);
    end

    dvp_Arq(dvp_Arq>1)=1;
    dvp_Harq(dvp_Harq>1)=1;

    figure(1)
    plot(1e3*delayTargetArray/rtt,log10(dvp_Arq),lineOrder(3),...
        'Marker',markerOrder(1),'MarkerSize',10,...
        'color',corder(rttIter,:),'LineWidth', 3);
    hold on;
    plot(1e3*delayTargetArray/rtt,log10(dvp_Harq),lineOrder(1),...
        'Marker',markerOrder(1),'MarkerSize',10,...
        'color',corder(rttIter,:),'LineWidth', 3);
end

figure(1)
set(gca,'Fontsize',28)
% titletext = "DVP(total) at $d = "+d*1000+"$ ms vs. allocated resources per slot.";
% title(titletext,'fontsize',36,'interpreter','latex');
grid on;
legend( 'ARQ: RTT$ = 1$ ms' ,'HARQ: RTT$ = 1$ ms' ,...
        'ARQ: RTT$ = 3$ ms' ,'HARQ: RTT$ = 3$ ms' ,...
        'ARQ: RTT$ = 5$ ms' ,'HARQ: RTT$ = 5$ ms' ,...
        'ARQ: RTT$ = 7$ ms' ,'HARQ: RTT$ = 7$ ms' ,...
    'fontsize',32,'interpreter','latex','location','northeast','NumColumns',2,'Orientation','horizontal');
xlabel('$d/\mathrm{RTT}$','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');
% grid minor
ylim([-10,0])
xlim([0,10])


figure(1)
saveas(gcf,'Figures/bler_vs_dByRttRatio.fig')
saveas(gcf,'Figures/bler_vs_dByRttRatio.epsc')
movefile('Figures/bler_vs_dByRttRatio.epsc','Figures/bler_vs_dByRttRatio.eps')




% HARQ saturates at ~1e-7 because of the transmission limit of 4.
% The RTT=1 plot behaves so because we considered threshold values which
% are 1ms apart and there was only one sample between each point. A better
% plot is given below
