


clear;clc;close all;
addpath('saved data')
addpath('funcs')

load 'All_Pop.mat'

R = All_Pop;
secretKey = 'dbSNP';



tau = 0.75;
Omega = 10^-5;
beta1 = 10^-5;
c0_list = [2:1:6];
c1_list = c0_list;

N = 15000;
T = 156;

epsilon_list = [];
Acc_gen_scope_tardos = [];
for i = 1:length(c0_list)
    c0  = c0_list(i);
    c1  = c1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);
    epsilon_list = [epsilon_list epsilon];
    [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1);
    acc_fixed_Omega = 0;
    for i = 1:c1
        MR = LDP_FP_dbSNP_one_seed(R,secretKey,epsilon,Tardos_code(i,:));
        MR_list{i} = MR;
        acc_fixed_Omega = acc_fixed_Omega + ...
            sum(sum(  MR(:,[2:end]).Variables==R(:,[2:end]).Variables  ))/ (N*T);
    end
    Acc_gen_scope_tardos = [Acc_gen_scope_tardos  acc_fixed_Omega/c1];
end

% save('varyOmega01.mat')

%%
Acc_Two_Stage = [];
for i = 1:length(c0_list)
    c0  = c0_list(i);
    c1  = c1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);

    acc_two_stage = 0;
    for i = 1:c1
        LDP_R = LDP_dbSNP(R,epsilon);
        [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1);
        MR = FP_dbSNP_one_seed(LDP_R,secretKey,gamma,Tardos_code);
        acc_two_stage = acc_two_stage + ...
            sum(sum(  MR(:,[2:end]).Variables==R(:,[2:end]).Variables  ))/ (N*T);
    end
    Acc_Two_Stage = [Acc_Two_Stage  acc_two_stage/c1];
end




