function [expectedBlerVector] = getLenaHarqBler_old( ...
        iterLen, mcsIndex, snrRcvd_avg, hvar, M_Harq, ...
        bgType, cbsTarget, Nrb)

    switch mcsIndex
        case 20
            expectedBlerVector = [0.874926, 0.014023, 0.005965, 0.0];
        case 14
            expectedBlerVector = [0.096412, 0.002761, 0.016807, 0.0];
        case 12
            expectedBlerVector = [0.041926, 0.003761, 0.000000, 0.0];
        case 9
            expectedBlerVector = [0.011856, 0.001126, 0.0, 0.0];
        case 8
            expectedBlerVector = [0.006413, 0.000637, 0.0, 0.0];
        case 7
            expectedBlerVector = [0.003211, 0.000317, 0.0, 0.0];
        case 6
            expectedBlerVector = [0.001405, 0.000724, 0.0, 0.0];
        case 5
            expectedBlerVector = [0.000658, 0.0, 0.0, 0.0];
        case 4
            expectedBlerVector = [0.000228, 0.0, 0.0, 0.0];
        case 3
            expectedBlerVector = [0.000084, 0.0, 0.0, 0.0];
        case 2
            expectedBlerVector = [0.000011, 0.0, 0.0, 0.0];
        case 1
            expectedBlerVector = [0.000011, 0.0, 0.0, 0.0];
        case 0
            expectedBlerVector = [0.000011, 0.0, 0.0, 0.0];
        otherwise
            error('No LENA table is available for MCS %d.', mcsIndex);
    end
% getLenaHarqBler
%   (MCS, mean SNR) -> HARQ 1:M_Harq BLER (5G LENA tables and EESM)
%
% Inputs:
%   iterLen     : Number of Monte Carlo iterations
%   mcsIndex    : 0~28
%   snrRcvd_avg : Mean received SNR (linear scale)
%   hvar        : Mean fading power (exprnd parameter)
%   M_Harq      : Maximum number of HARQ attempts
%   bgType      : 1 or 2 (BG TYPE)
%   cbsTarget   : Target CBS (e.g., 816; selects the nearest available CBS)
%   betaEesm    : EESM beta parameter
%
% Outputs:
%   expectedBlerVector : 1 x M_Harq, mean BLER after each HARQ attempt

    % blerArray   = zeros(iterLen, M_Harq);
    % fadingCoeff = normrnd(0,2.6,M_Harq,iterLen);%exprnd(hvar, M_Harq, iterLen);  % [M_Harq x iterLen]    
    % % beta = getBetaTable1(mcsIndex);
    % 
    % 
    % for iter = 1:iterLen
    %     snrRcvd_inst = (snrRcvd_avg)+fadingCoeff(:,iter);
    %     for m = 1:M_Harq
    %         sinrEff_dB = 10*log10(mean(snrRcvd_inst(1:m))*m);
    % 
    %         % sinrEff_lin = -beta*log(mean(exp(-gamma/beta)));          % Sum in the linear domain
    %         % % sinrEff_lin
    %         % sinrEff_lin = max(sinrEff_lin, 1e-12);
    %         % sinrEff_dB  = 10*log10(sinrEff_lin)+10*log10(Nrb);
    %         % 
    %         % ---- LENA table-based BLER ----
    %         bler = getLenaBlerFromSinr(mcsIndex, sinrEff_dB, bgType);
    % 
    %         blerArray(iter, m) = bler;
    %     end
    % end
    % 
    % expectedBlerVector = mean(blerArray, 1);  % [1 x M_Harq]
