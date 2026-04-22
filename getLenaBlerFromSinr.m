function bler = getLenaBlerFromSinr(mcsIndex, sinr_dB, bgType)
% getLenaBlerFromSinr
%   nr-eesm-t1.cc에서 생성한 테이블 기반으로
%   (MCS, SINR[dB]) -> BLER 매핑
%
% 입력:
%   mcsIndex : 0 ~ 28
%   sinr_dB  : 스칼라 또는 벡터 (dB 단위)
%   bgType   : 1 or 2 (BG TYPE 1 / 2)
%   cbsTarget: 아무 값이나 넣어도 됨 (가장 가까운 CBS 자동 선택)
%
% 출력:
%   bler     : sinr_dB와 같은 크기의 BLER 값

    % 1) BG 타입에 따라 테이블 함수 선택
    switch bgType
        case 1
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG1(mcsIndex);
        case 2
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG2(mcsIndex);
        otherwise
            error('bgType must be 1 or 2');
    end

    % 2) 보간으로 BLER 계산 (벡터 SINR도 가능)
    % log-domain 보간이 더 안정적이라면 이렇게 해도 됨:
    %   logBler = interp1(..., log10(blerTable + 1e-12), ...);
    %   bler    = 10.^logBler;
    bler = interp1( ...
        snrTable_dB, blerTable, ...
        sinr_dB, 'linear', 'extrap');

    % 0~1 사이로 clip
    bler = min(max(bler, 0.0), 1.0);
end
