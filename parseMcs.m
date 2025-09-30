function [Qm, codeRate, specEff] = parseMcs(mcs,mcsTableNum)
if(isempty(mcs))
    error("ERROR: MCS is a mandatory input")
end
if nargin < 2
    mcsTableNum = 1;
end

if mcsTableNum == 1
    % 38.214  - Table 5.1.3.1-1: MCS index table 1 for PDSCH
    switch mcs
        case 0
            Qm = 2; codeRate = 120; specEff = 0.2344;
        case 1
            Qm = 2; codeRate = 157; specEff = 0.3066;
        case 2
            Qm = 2; codeRate = 193; specEff = 0.3770;
        case 3
            Qm = 2; codeRate = 251; specEff = 0.4902;
        case 4
            Qm = 2; codeRate = 308; specEff = 0.6016;
        case 5
            Qm = 2; codeRate = 379; specEff = 0.7402;
        case 6
            Qm = 2; codeRate = 449; specEff = 0.8770;
        case 7
            Qm = 2; codeRate = 526; specEff = 1.0273;
        case 8
            Qm = 2; codeRate = 602; specEff = 1.1758;
        case 9
            Qm = 2; codeRate = 679; specEff = 1.3262;
        case 10
            Qm = 4; codeRate = 340; specEff = 1.3281;
        case 11
            Qm = 4; codeRate = 378; specEff = 1.4766;
        case 12
            Qm = 4; codeRate = 434; specEff = 1.6953;
        case 13
            Qm = 4; codeRate = 490; specEff = 1.9141;
        case 14
            Qm = 4; codeRate = 553; specEff = 2.1602;
        case 15
            Qm = 4; codeRate = 616; specEff = 2.4063;
        case 16
            Qm = 4; codeRate = 658; specEff = 2.5703;
        case 17
            Qm = 6; codeRate = 438; specEff = 2.5664;
        case 18
            Qm = 6; codeRate = 466; specEff = 2.7305;
        case 19
            Qm = 6; codeRate = 517; specEff = 3.0293;
        case 20
            Qm = 6; codeRate = 567; specEff = 3.3223;
        case 21
            Qm = 6; codeRate = 616; specEff = 3.6094;
        case 22
            Qm = 6; codeRate = 666; specEff = 3.9023;
        case 23
            Qm = 6; codeRate = 719; specEff = 4.2129;
        case 24
            Qm = 6; codeRate = 772; specEff = 4.5234;
        case 25
            Qm = 6; codeRate = 822; specEff = 4.8164;
        case 26
            Qm = 6; codeRate = 873; specEff = 5.1152;
        case 27
            Qm = 6; codeRate = 910; specEff = 5.3320;
        case 28
            Qm = 6; codeRate = 948; specEff = 5.5547;
        otherwise
            Qm = -1; codeRate = -1; specEff = -1;
    end
elseif mcsTableNum == 2
    % 38.214 - Table 5.1.3.1-2: MCS index table 2 for PDSCH
    switch mcs
        case 0
            Qm = 2;  codeRate = 120  ;  specEff = 0.2344;
        case 1
            Qm = 2;  codeRate = 193  ;  specEff = 0.377;
        case 2
            Qm = 2;  codeRate = 308  ;  specEff = 0.6016;
        case 3
            Qm = 2;  codeRate = 449  ;  specEff = 0.877;
        case 4
            Qm = 2;  codeRate = 602  ;  specEff = 1.1758;
        case 5
            Qm = 4;  codeRate = 378  ;  specEff = 1.4766;
        case 6
            Qm = 4;  codeRate = 434  ;  specEff = 1.6953;
        case 7
            Qm = 4;  codeRate = 490  ;  specEff = 1.9141;
        case 8
            Qm = 4;  codeRate = 553  ;  specEff = 2.1602;
        case 9
            Qm = 4;  codeRate = 616  ;  specEff = 2.4063;
        case 10
            Qm = 4;  codeRate = 658  ;  specEff = 2.5703;
        case 11
            Qm = 6;  codeRate = 466  ;  specEff = 2.7305;
        case 12
            Qm = 6;  codeRate = 517  ;  specEff = 3.0293;
        case 13
            Qm = 6;  codeRate = 567  ;  specEff = 3.3223;
        case 14
            Qm = 6;  codeRate = 616  ;  specEff = 3.6094;
        case 15
            Qm = 6;  codeRate = 666  ;  specEff = 3.9023;
        case 16
            Qm = 6;  codeRate = 719  ;  specEff = 4.2129;
        case 17
            Qm = 6;  codeRate = 772  ;  specEff = 4.5234;
        case 18
            Qm = 6;  codeRate = 822  ;  specEff = 4.8164;
        case 19
            Qm = 6;  codeRate = 873  ;  specEff = 5.1152;
        case 20
            Qm = 8;  codeRate = 682.5;  specEff = 5.332;
        case 21
            Qm = 8;  codeRate = 711  ;  specEff = 5.5547;
        case 22
            Qm = 8;  codeRate = 754  ;  specEff = 5.8906;
        case 23
            Qm = 8;  codeRate = 797  ;  specEff = 6.2266;
        case 24
            Qm = 8;  codeRate = 841  ;  specEff = 6.5703;
        case 25
            Qm = 8;  codeRate = 885  ;  specEff = 6.9141;
        case 26
            Qm = 8;  codeRate = 916.5;  specEff = 7.1602;
        case 27
            Qm = 8;  codeRate = 948  ;  specEff = 7.4063;
        otherwise
            Qm = -1; codeRate = -1; specEff = -1;
    end
elseif mcsTableNum == 3
    % 38.214 - Table 5.1.3.1-3: MCS index table 3 for PDSCH
    switch mcs
        case 0
            Qm = 2;   codeRate = 30 ;   specEff = 0.0586;
        case 1
            Qm = 2;   codeRate = 40 ;   specEff = 0.0781;
        case 2
            Qm = 2;   codeRate = 50 ;   specEff = 0.0977;
        case 3
            Qm = 2;   codeRate = 64 ;   specEff = 0.125;
        case 4
            Qm = 2;   codeRate = 78 ;   specEff = 0.1523;
        case 5
            Qm = 2;   codeRate = 99 ;   specEff = 0.1934;
        case 6
            Qm = 2;   codeRate = 120;   specEff = 0.2344;
        case 7
            Qm = 2;   codeRate = 157;   specEff = 0.3066;
        case 8
            Qm = 2;   codeRate = 193;   specEff = 0.377;
        case 9
            Qm = 2;   codeRate = 251;   specEff = 0.4902;
        case 10
            Qm = 2;   codeRate = 308;   specEff = 0.6016;
        case 11
            Qm = 2;   codeRate = 379;   specEff = 0.7402;
        case 12
            Qm = 2;   codeRate = 449;   specEff = 0.877;
        case 13
            Qm = 2;   codeRate = 526;   specEff = 1.0273;
        case 14
            Qm = 2;   codeRate = 602;   specEff = 1.1758;
        case 15
            Qm = 4;   codeRate = 340;   specEff = 1.3281;
        case 16
            Qm = 4;   codeRate = 378;   specEff = 1.4766;
        case 17
            Qm = 4;   codeRate = 434;   specEff = 1.6953;
        case 18
            Qm = 4;   codeRate = 490;   specEff = 1.9141;
        case 19
            Qm = 4;   codeRate = 553;   specEff = 2.1602;
        case 20
            Qm = 4;   codeRate = 616;   specEff = 2.4063;
        case 21
            Qm = 6;   codeRate = 438;   specEff = 2.5664;
        case 22
            Qm = 6;   codeRate = 466;   specEff = 2.7305;
        case 23
            Qm = 6;   codeRate = 517;   specEff = 3.0293;
        case 24
            Qm = 6;   codeRate = 567;   specEff = 3.3223;
        case 25
            Qm = 6;   codeRate = 616;   specEff = 3.6094;
        case 26
            Qm = 6;   codeRate = 666;   specEff = 3.9023;
        case 27
            Qm = 6;   codeRate = 719;   specEff = 4.2129;
        case 28
            Qm = 6;   codeRate = 772;   specEff = 4.5234;
        otherwise
            Qm = -1; codeRate = -1; specEff = -1;
    end
else
    error("ERROR: Unsupported MCS Table number.");
end
codeRate = codeRate/1024;
end

