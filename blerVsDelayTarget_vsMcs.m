clearvars -except h1 h2 h3;clc

corder = colororder;
lineOrder = ["-","--",":","-","--",":","-."];
markerOrder = ["none","o","x","diamond","|","square"];

isReTxPrioritised = true;
isSimRequired = false;
isArqRequired = false;

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
f=1/3;
delta_Orig = 2;
pktLen = 8 * 100; % Packet length in bits


delayTargetArrayOrig  = (2.5:0.1:25) * 1e-3;
McsArray = 10:2:20;

zeta = 1; % Decoding delay. A parameter designed later in the research.
delta = delta_Orig+zeta;
delayTargetArray=delayTargetArrayOrig-zeta*1e-3;
if(min(delayTargetArray)<1e-3)
    error("Invalid delayTarget,zeta combination. Correct it!")
end


delayTargetLen = length(delayTargetArray);
for mccIdx = 1:length(McsArray)
    mcs = McsArray(mccIdx)

    [~, ~, eta] = parseMcs(mcs,1);
    Nrb = ceil(pktLen/(180*eta));
    Nre = Nrb *180;
    B = 180e3 * Nrb;
    dataBits = Nre *eta;

    dvp_Arq  = zeros(1,delayTargetLen);
    dvp_Arq_sim = zeros(1,delayTargetLen);
    dvp_Harq  = zeros(1,delayTargetLen);
    dvp_Harq_sim = zeros(1,delayTargetLen);

    for dIter = 1:delayTargetLen
        d = delayTargetArray(dIter);
        if(isArqRequired)


            pArq = getFblPer(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar);


            kd = floor(((d/T) + delta)/(1 + delta));

            q = pArq/(1 - f);

            term1 = pArq^kd;
            term2 = f * (q)^(floor(d/T) + delta);
            term3 = 1 - pArq^kd * (q)^(-1 * kd * (1 + delta));
            term4 = 1 - f - (q)^delta;

            dvp_Arq(dIter) = term1 - (term2 * term3) / term4;


            if(isSimRequired)
                [~,~,latencyInSlotsVector,~] = simulate_queue_noIR(isReTxPrioritised,f,pArq(nrbIter),delta,simLen,Qmax,M,false);
                delaytargetInSlots = (d/T);
                dvp_Arq_sim(dIter) = (sum(latencyInSlotsVector > delaytargetInSlots))/length(latencyInSlotsVector);
            end
        end

        p_Harq = getFblPerVector_Harq(simLen_MonteCarlo,snr_avgRcvd,eta,Nre,V,hvar,M_Harq);

        dvp_Harq(dIter) = getDvpHarq (d, delta, p_Harq, f, T, M_Harq, Qmax_Harq);

    end

    dvp_Arq(dvp_Arq>1)=1;
    dvp_Harq(dvp_Harq>1)=1;

    if(isArqRequired)

        if(isSimRequired)
            plot(delayTargetArrayOrig*1000,log10(dvp_Arq),lineOrder(mccIdx),'Marker',markerOrder(mccIdx),'MarkerSize',10,'color',corder(1,:),'LineWidth', 3);
            hold on
            plot(delayTargetArrayOrig*1000,log10(dvp_Arq_sim),lineOrder(mccIdx),'Marker',markerOrder(mccIdx),'MarkerSize',10,'color',corder(2,:),'LineWidth', 3);

        else
            plot(delayTargetArrayOrig*1000,log10(dvp_Arq),lineOrder(3),...
                'Marker',markerOrder(mccIdx),'MarkerSize',2,...
                'color',corder(mccIdx,:),'LineWidth', 3);
            hold on;
            plot(delayTargetArrayOrig*1000,log10(dvp_Harq),lineOrder(1),...
                'Marker',markerOrder(mccIdx),'MarkerSize',2,...
                'color',corder(mccIdx,:),'LineWidth', 3);
        end
    else
        plot(delayTargetArrayOrig*1000,log10(dvp_Harq),lineOrder(1),...
            'Marker',markerOrder(mccIdx),'MarkerSize',2,...
            'color',corder(mccIdx,:),'LineWidth', 3); hold on
    end
    drawnow
end

set(gca,'Fontsize',28)
grid on;
if(isArqRequired)
    legend( 'ARQ: MCS$ = 12$' ,'HARQ: MCS$ = 12$' ,...
        'ARQ: MCS$ = 14$' ,'HARQ: MCS$ = 14$' ,...
        'ARQ: MCS$ = 16$' ,'HARQ: MCS$ = 16$' ,...
        'ARQ: MCS$ = 18$' ,'HARQ: MCS$ = 18$' ,...
        'ARQ: MCS$ = 20$' ,'HARQ: MCS$ = 20$',...
        'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');
else
    legend( 'HARQ: MCS$ = 12$' ,...
        'HARQ: MCS$ = 14$' ,...
        'HARQ: MCS$ = 16$' ,...
        'HARQ: MCS$ = 18$' ,...
        'HARQ: MCS$ = 20$',...
        'fontsize',32,'interpreter','latex','location','southwest','NumColumns',2,'Orientation','horizontal');

end
xlabel('Delay Target (ms)','fontsize',42,'interpreter','latex');
ylabel('log$_{10}$(DVP)','fontsize',42,'interpreter','latex');


title('DVP vs. $d$: $f=1/3, \delta=2,\zeta=1, n=800, T=1ms$','fontsize',42,'interpreter','latex')