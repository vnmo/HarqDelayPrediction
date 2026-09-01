function bler = getLenaBlerFromSinr(mcsIndex, sinr_dB, bgType)

    switch bgType
        case 1
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG1(mcsIndex);
        case 2
            [snrTable_dB, blerTable, ~] = getLenaTable_EESM_BG2(mcsIndex);
        otherwise
            error('bgType must be 1 or 2');
    end

    bler = interp1( ...
        snrTable_dB, blerTable, ...
        sinr_dB, 'linear', 'extrap');

    bler = min(max(bler, 0.0), 1.0);
end
