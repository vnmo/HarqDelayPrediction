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

