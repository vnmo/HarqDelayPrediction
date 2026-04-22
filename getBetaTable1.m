function beta = getBetaTable1(mcs)
% getBetaTable1 : return optimal beta (Table 1) for given MCS (0~28)

switch mcs
    case 0
        beta = 1.60;
    case 1
        beta = 1.61;
    case 2
        beta = 1.63;
    case 3
        beta = 1.65;
    case 4
        beta = 1.67;
    case 5
        beta = 1.70;
    case 6
        beta = 1.73;
    case 7
        beta = 1.76;
    case 8
        beta = 1.79;
    case 9
        beta = 1.82;
    case 10
        beta = 3.97;
    case 11
        beta = 4.27;
    case 12
        beta = 4.71;
    case 13
        beta = 5.16;
    case 14
        beta = 5.66;
    case 15
        beta = 6.16;
    case 16
        beta = 6.50;
    case 17
        beta = 9.95;
    case 18
        beta = 10.97;
    case 19
        beta = 12.92;
    case 20
        beta = 14.96;
    case 21
        beta = 17.06;
    case 22
        beta = 19.33;
    case 23
        beta = 21.85;
    case 24
        beta = 24.51;
    case 25
        beta = 27.14;
    case 26
        beta = 29.94;
    case 27
        beta = 32.05;
    case 28
        beta = 34.28;
    otherwise
        error('MCS must be an integer between 0 and 28.');
end
end
