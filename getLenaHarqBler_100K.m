function [expectedBlerVector] = getLenaHarqBler( ...
        iterLen, mcsIndex, snrRcvd_avg, hvar, M_Harq, ...
        bgType, cbsTarget, Nrb)

sinrDb = round(snrRcvd_avg); % assume caller passes ~folder SINR(dB)
expectedBlerVector = [NaN, NaN, NaN, NaN];

switch sinrDb
    case 10
        switch mcsIndex
            case 1
                expectedBlerVector = [0.000202, 0.536762, 0.889196, 0.512668];
            case 2
                expectedBlerVector = [0.000584, 0.263832, 0.897722, 0.548796];
            case 3
                expectedBlerVector = [0.001440, 0.259153, 0.747827, 0.502124];
            case 4
                expectedBlerVector = [0.002679, 0.253640, 0.551185, 0.497803];
            case 5
                expectedBlerVector = [0.004430, 0.179492, 0.678575, 0.448440];
            case 6
                expectedBlerVector = [0.006402, 0.211142, 0.456250, 0.446178];
            case 7
                expectedBlerVector = [0.009085, 0.232585, 0.500156, 0.396729];
            case 8
                expectedBlerVector = [0.012088, 0.186822, 0.708696, 0.322709];
            case 9
                expectedBlerVector = [0.013943, 0.219580, 0.575499, 0.453297];
            case 12
                expectedBlerVector = [0.018922, 0.479023, 0.923590, 0.782191];
            case 14
                expectedBlerVector = [0.040827, 0.229478, 0.920352, 0.748815];
            case 20
                expectedBlerVector = [0.333387, 0.235711, 0.867227, 0.848201];
            otherwise
                error('No table is available for SINR=%d dB, MCS %d.', sinrDb, mcsIndex);
        end

    case 20
        switch mcsIndex
            case 1
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 2
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 3
                expectedBlerVector = [0.000013, 0.000000, 0.000000, 0.000000];
            case 4
                expectedBlerVector = [0.000088, 0.031956, 0.000000, 0.000000];
            case 5
                expectedBlerVector = [0.000183, 0.020797, 0.003400, 0.000000];
            case 6
                expectedBlerVector = [0.000299, 0.156078, 0.000000, 0.000000];
            case 7
                expectedBlerVector = [0.000514, 0.162687, 0.125911, 0.275212];
            case 8
                expectedBlerVector = [0.000993, 0.101147, 0.475519, 0.080119];
            case 9
                expectedBlerVector = [0.001352, 0.097631, 0.567039, 0.229047];
            case 12
                expectedBlerVector = [0.001103, 0.438580, 0.931162, 0.659335];
            case 14
                expectedBlerVector = [0.002471, 0.204121, 0.935274, 0.744152];
            case 20
                expectedBlerVector = [0.008586, 0.288354, 0.908075, 0.811415];
            otherwise
                error('No table is available for SINR=%d dB, MCS %d.', sinrDb, mcsIndex);
        end

    case 30
        switch mcsIndex
            case 1
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 2
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 3
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 4
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 5
                expectedBlerVector = [0.000000, 0.000000, 0.000000, 0.000000];
            case 6
                expectedBlerVector = [0.000014, 0.000000, 0.000000, 0.000000];
            case 7
                expectedBlerVector = [0.000014, 0.000000, 0.000000, 0.000000];
            case 8
                expectedBlerVector = [0.000044, 0.000000, 0.000000, 0.000000];
            case 9
                expectedBlerVector = [0.000165, 0.000000, 0.000000, 0.000000];
            case 12
                expectedBlerVector = [0.000065, 0.136564, 0.811719, 0.270573];
            case 14
                expectedBlerVector = [0.000299, 0.034551, 1.000000, 0.334067];
            case 20
                expectedBlerVector = [0.000803, 0.205123, 0.961898, 0.777231];
            otherwise
                error('No table is available for SINR=%d dB, MCS %d.', sinrDb, mcsIndex);
        end

    otherwise
        error('No table is available for SINR=%d dB.', sinrDb);
end
end
