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
            error('MCS %d에 대한 LENA 테이블이 구현되지 않았습니다.', mcsIndex);
    end
% getLenaHarqBler
%   (MCS, 평균 SNR) -> HARQ 1~M_Harq BLER (5G LENA 테이블 + EESM 기반)
%
% 입력:
%   iterLen     : Monte Carlo 반복 수
%   mcsIndex    : 0~28
%   snrRcvd_avg : 평균 수신 SNR (선형 값)
%   hvar        : Rayleigh 평균 (exprnd 파라미터)
%   M_Harq      : HARQ 최대 전송 횟수
%   bgType      : 1 or 2 (BG TYPE)
%   cbsTarget   : CBS 타겟 (예: 816, 아무거나 넣어도 가장 가까운 CBS 선택)
%   betaEesm    : EESM beta 파라미터
%
% 출력:
%   expectedBlerVector : 1 x M_Harq, 각 HARQ r회까지의 평균 BLER

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
    %         % sinrEff_lin = -beta*log(mean(exp(-gamma/beta)));          % 선형 도메인 합
    %         % % sinrEff_lin
    %         % sinrEff_lin = max(sinrEff_lin, 1e-12);
    %         % sinrEff_dB  = 10*log10(sinrEff_lin)+10*log10(Nrb);
    %         % 
    %         % ---- LENA 테이블 기반 BLER ----
    %         bler = getLenaBlerFromSinr(mcsIndex, sinrEff_dB, bgType);
    % 
    %         blerArray(iter, m) = bler;
    %     end
    % end
    % 
    % expectedBlerVector = mean(blerArray, 1);  % [1 x M_Harq]