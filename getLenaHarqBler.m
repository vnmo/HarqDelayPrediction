function [expectedBlerVector] = getLenaHarqBler( ...
        iterLen, mcsIndex, snrRcvd_avg, hvar, M_Harq, ...
        bgType, cbsTarget, Nrb)

sinrDb = round(snrRcvd_avg); % assume caller passes ~folder SINR(dB)
expectedBlerVector = [NaN, NaN, NaN, NaN];

switch sinrDb
    case 10
        switch mcsIndex
            case 0
                expectedBlerVector = [0.000129, 0.805690, 0.865164, 0.473989];
            case 1
                expectedBlerVector = [0.000292, 0.512591, 0.893291, 0.515929];
            case 2
                expectedBlerVector = [0.000657, 0.336892, 0.897005, 0.523517];
            case 3
                expectedBlerVector = [0.001541, 0.269473, 0.738805, 0.495690];
            case 4
                expectedBlerVector = [0.002802, 0.255357, 0.538072, 0.511996];
            case 5
                expectedBlerVector = [0.004669, 0.183274, 0.708992, 0.432244];
            case 6
                expectedBlerVector = [0.006788, 0.213701, 0.462340, 0.445359];
            case 7
                expectedBlerVector = [0.009277, 0.238621, 0.477600, 0.424872];
            case 8
                expectedBlerVector = [0.011856, 0.201959, 0.713967, 0.304975];
            case 9
                expectedBlerVector = [0.013527, 0.237267, 0.577270, 0.452832];
            case 12
                expectedBlerVector = [0.018851, 0.467728, 0.923449, 0.793764];
            case 14
                expectedBlerVector = [0.041093, 0.222699, 0.921609, 0.762522];
            case 20
                expectedBlerVector = [0.338911, 0.233558, 0.867235, 0.849936];
            otherwise
                error('Table Error [No Table]: SINR=%d dB, MCS %d', sinrDb, mcsIndex);
        end
    case 20
        switch mcsIndex
            case 0
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 1
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 2
                expectedBlerVector = [0.000001, 0.000000, 0.000000, 0.000000];
            case 3
                expectedBlerVector = [0.000002, 0.306710, 0.684127, 0.399897];
            case 4
                expectedBlerVector = [0.000008, 0.130160, 0.694898, 0.427881];
            case 5
                expectedBlerVector = [0.000020, 0.094430, 0.492896, 0.436133];
            case 6
                expectedBlerVector = [0.000034, 0.138910, 0.221728, 0.431658];
            case 7
                expectedBlerVector = [0.000069, 0.124272, 0.371064, 0.294501];
            case 8
                expectedBlerVector = [0.000125, 0.078968, 0.634785, 0.205617];
            case 9
                expectedBlerVector = [0.000166, 0.105067, 0.413849, 0.341545];
            case 12
                expectedBlerVector = [0.000165, 0.386781, 0.894452, 0.700177];
            case 14
                expectedBlerVector = [0.000369, 0.181887, 0.920521, 0.708970];
            case 20
                expectedBlerVector = [0.003121, 0.148150, 0.884651, 0.765015];
            otherwise
                error('Table Error [No Table]: SINR=%d dB, MCS %d', sinrDb, mcsIndex);
        end

    

    otherwise
        error('Table Error [No Table]: SINR=%d dB', sinrDb);
end
end
