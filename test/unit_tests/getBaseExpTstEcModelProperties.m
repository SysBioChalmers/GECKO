function [expRxns, expMetNames, expS] = getBaseExpTstEcModelProperties()
    expRxns = {'R1';'R1_REV';'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_2';'R3';'R4';'R5';'S1';'S2'; ...
        'usage_prot_P1';'usage_prot_P2';'usage_prot_P3';'usage_prot_P4';'usage_prot_P5';'prot_pool_exchange'};
    expMetNames = {'e1';'e2';'m1';'m2';'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5';'prot_pool'};

    expS = sparse(zeros(length(expMetNames), length(expRxns)));
    expS(1,1) = -1;%R1
    expS(3,1) = 1;
    expS(1,2) = 1;%R1_REV
    expS(3,2) = -1;
    expS(3,3) = -1;%R2_EXP_1
    expS(4,3) = 1;
    expS(3,4) = -1;%R2_EXP_2
    expS(4,4) = 1;
    expS(3,5) = 1;%R2_REV_EXP_1
    expS(4,5) = -1;
    expS(3,6) = 1;%R2_REV_EXP_2
    expS(4,6) = -1;
    expS(3,7) = -1;%R3
    expS(4,7) = 1;
    expS(3,8) = -1;%R4
    expS(4,8) = 1;
    expS(4,9) = -1;%R5
    expS(2,9) = 1;
    expS(1,10) = -1;%E1
    expS(2,11) = -1;%E2
    %Now add the protein reactions
    expS(10,12) = -1;%P1
    expS(5,12) = 1;
    expS(10,13) = -1;%P2
    expS(6,13) = 1;
    expS(10,14) = -1;%P3
    expS(7,14) = 1;
    expS(10,15) = -1;%P4
    expS(8,15) = 1;
    expS(10,16) = -1;%P5
    expS(9,16) = 1;
    expS(10,17) = 1;%prot pool
end
