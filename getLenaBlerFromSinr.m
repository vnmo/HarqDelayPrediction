function bler = getLenaBlerFromSinr(mcsIndex, sinr_dB, bgType)
% getLenaBlerFromSinr
%   Based on tables generated from nr-eesm-t1.cc
%   (MCS, SINR[dB]) -> BLER mapping
%
% Inputs:
%   mcsIndex : 0 ~ 28
%   sinr_dB  : Scalar or vector in dB
%   bgType   : 1 or 2 (BG TYPE 1 / 2)
%   CBS selection is handled by the table function.
%
% Outputs:
%   bler     : BLER values with the same shape as sinr_dB

    % 1) Select the table function by BG type
    switch bgType
        case 1
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG1(mcsIndex);
        case 2
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG2(mcsIndex);
        otherwise
            error('bgType must be 1 or 2');
    end

    % 2) Interpolate BLER (also supports vector SINR)
    % Alternative log-domain interpolation:
    %   logBler = interp1(..., log10(blerTable + 1e-12), ...);
    %   bler    = 10.^logBler;
    bler = interp1( ...
        snrTable_dB, blerTable, ...
        sinr_dB, 'linear', 'extrap');

    % Clip to [0,1]
    bler = min(max(bler, 0.0), 1.0);
end
