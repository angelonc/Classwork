%% PROBABILISTIC MODELING HW3

%% B
% B.1
fprintf('\nB.1\n\n')
% Conditions [A B C]:
conditions = [0 0 0; ...
              0 0 1; ...
              0 1 0; ...
              0 1 1; ...
              1 0 0; ...
              1 0 1; ...
              1 1 0; ...
              1 1 1];

joints =     [.192; ...
              .144; ...
              .048; ...
              .216; ...
              .192; ...
              .064; ...
              .048; ...
              .096];

conditions = logical(conditions);

% Find p(A,B) by summing over C:
pAB = sum(joints(conditions(:,1) & conditions(:,2)));

% Find p(A), p(B) by summing over other variables, respectively:
pA = sum(joints(conditions(:,1)));
pB = sum(joints(conditions(:,2)));

% Dependent if p(A,B) != p(A)*p(B);
tmp = pAB == pA * pB;
fprintf('Does p(A,B) == p(A)p(B)? (0 = no, 1 = yes): %g\n',tmp);
if tmp
    fprintf('\tIndependent\n');
else
    fprintf('\tDependent\n');
end

% Condition over C for C == 0 and 1
for i = [0 1]
    pC = sum(joints(conditions(:,3) == i));
    pABC = joints(conditions(:,1) & conditions(:,2) ...
                   & conditions(:,3) == i);
    
    pAB_C = pABC / pC;
    
    pAC = sum(joints(conditions(:,1) & conditions(:,3) == i));
    pBC = sum(joints(conditions(:,2) & conditions(:,3) == i));
    
    pA_C = pAC / pC;
    pB_C = pBC / pC;
    
    tmp = pAB_C == pA_C * pB_C;
    fprintf('Does p(A,B|C) == p(A|C)p(B|C) when C is %g?: %g\n',i,tmp);
    if tmp
        fprintf('\tIndependent\n');
    else
        fprintf('\tDependent\n');
    end
end

% B.2
fprintf('\nB.2\n\n');
pC_A = pAC / pA;
pABC = joints(conditions(:,1) & conditions(:,2) ...
                   & conditions(:,3));

tmp = pABC == pA * pC_A * pB_C;
fprintf('Does p(A,B,C) == p(A)p(C|A)p(B|C)?: %g\n\n',tmp);
clear all

%% C
% Probabilities
pB1      = .9;
pB0      = .1;
pF1      = .9;
pF0      = .1;

pG1_B1F1 = .8;
pG1_B1F0 = .2;
pG1_B0F1 = .2;
pG1_B0F0 = .1;
pG0_B1F1 = .2;
pG0_B1F0 = .8;
pG0_B0F1 = .8;
pG0_B0F0 = .9;

pD1_G1   = .9;
pD0_G1   = .1;
pD0_G0   = .9;
pD1_G0   = .1;

% Conditions [B F G D]:
conditions = [0 0 0 0; ...
              0 0 1 0; ...
              0 1 0 0; ...
              0 1 1 0; ...
              1 0 0 0; ...
              1 0 1 0; ...
              1 1 0 0; ...
              1 1 1 0];

% Calculate joints (p(B)p(F)p(G|B,F)p(D|G))... ugh:
joints(1)  = pB0*pF0*pG0_B0F0*pD0_G0;
joints(2)  = pB0*pF0*pG1_B0F0*pD0_G1;
joints(3)  = pB0*pF1*pG0_B0F1*pD0_G0;
joints(4)  = pB0*pF1*pG1_B0F1*pD0_G1;
joints(5)  = pB1*pF0*pG0_B1F0*pD0_G0;
joints(6)  = pB1*pF0*pG1_B1F0*pD0_G1;
joints(7)  = pB1*pF1*pG0_B1F1*pD0_G0;
joints(8)  = pB1*pF1*pG1_B1F1*pD0_G1;
% And the rest for posterities sake... ugh:
joints(9)  = pB0*pF0*pG0_B0F0*pD1_G0;
joints(10) = pB0*pF0*pG1_B0F0*pD1_G1;
joints(11) = pB0*pF1*pG0_B0F1*pD1_G0;
joints(12) = pB0*pF1*pG1_B0F1*pD1_G1;
joints(13) = pB1*pF0*pG0_B1F0*pD1_G0;
joints(14) = pB1*pF0*pG1_B1F0*pD1_G1;
joints(15) = pB1*pF1*pG0_B1F1*pD1_G0;
joints(16) = pB1*pF1*pG1_B1F1*pD1_G1;

% Sum to 1?
if sum(joints) ~= 1 error('Joints ~= 1'); end
joints = joints(1:8);

% C.1
pF0D0 = sum(joints(conditions(:,2) == 0 & conditions(:,4) == 0));
pD0 = sum(joints(conditions(:,4) == 0));
pF0_D0 = pF0D0 / pD0;


% C.2
pF0D0B0 = sum(joints(conditions(:,1) == 0 & conditions(:,2) == 0 & ...
                     conditions(:,4) == 0));
pD0B0 = sum(joints(conditions(:,1) == 0 & conditions(:,4) == 0));
pF0_D0B0 = pF0D0B0 / pD0B0;


