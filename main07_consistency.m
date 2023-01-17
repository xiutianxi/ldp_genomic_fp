


clear;clc;close all;
addpath('saved data')
addpath('funcs')

load 'All_Pop.mat'

R = All_Pop;

rng(10)

N = 15000;

T = 156;

temp = rand(N,1);

control_index = find(temp<0.5);
case_index = find(temp>0.5);

p_values = snp_p_value(R,control_index, case_index);

secretKey = 'dbSNP';

sp_id = 10;



Omega_list = 10.^[-13:1:-4];

tau = 0.85;
beta1_list  = 10.^[-13:1:-4];
c0 = 5;
c1 = c0;

topK = [10 20 30]/100;

[~,ground_truth_idx] = sort(p_values,'descend');

epsilon_list = [];


%% gen_scope
CONSIS_gen_scope = [];
for i = 1:length(Omega_list)
    Omega = Omega_list(i);
    beta1 = beta1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);
    epsilon_list = [epsilon_list epsilon];
    [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1);
    CONSIS = [];
    for t = 1:3
        K = ceil( topK(t) * 156 );
        consistency = 0;
        for j = 1:c1
            MR = LDP_FP_dbSNP_one_seed(R,secretKey,epsilon,Tardos_code(j,:));
            p_values = snp_p_value(MR,control_index, case_index);
            [~,idx] = sort(p_values,'descend');
            consistency = consistency + length( intersect (idx([1:K]),  ...
                ground_truth_idx([1:K])     )      )/K;
        end
        CONSIS = [CONSIS; consistency/c1];
    end
    CONSIS_gen_scope = [CONSIS_gen_scope CONSIS];
end

%% ldp
CONSIS_ldp = [];
for i = 1:length(Omega_list)
    Omega = Omega_list(i);
    beta1 = beta1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);
    CONSIS = [];
    for t = 1:3
        K = ceil( topK(t) * 156 );
        LDP_R = LDP_dbSNP(R,epsilon);
        p_values = snp_p_value(LDP_R,control_index, case_index);
        [~,idx] = sort(p_values,'descend');
        consistency =   length(  intersect ( idx([1:K]), ground_truth_idx([1:K])  )      )/K;
        CONSIS = [CONSIS; consistency];
    end
    CONSIS_ldp = [CONSIS_ldp CONSIS];
end
%% fp only
CONSIS_fp_only = [];
for i = 1:length(Omega_list)
    Omega = Omega_list(i);
    beta1 = beta1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);
    [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1);
    CONSIS = [];
    for t = 1:3
        K = ceil( topK(t) * 156 );
        consistency = 0;
        for j = 1:c1
            MR = FP_dbSNP_one_seed(R,secretKey,gamma,Tardos_code(j,:));
            p_values = snp_p_value(MR,control_index, case_index);
            [~,idx] = sort(p_values,'descend');
            consistency = consistency + length( intersect (idx([1:K]),  ...
                ground_truth_idx([1:K])     )      )/K;
        end
        CONSIS = [CONSIS; consistency/c1];
    end
    CONSIS_fp_only = [CONSIS_fp_only CONSIS];
end

%% two stage: ldp followed by fp
CONSIS_two_stage = [];
for i = 1:length(Omega_list)
    Omega = Omega_list(i);
    beta1 = beta1_list(i);
    gamma = 4 * pi^2 * c0^2 * log(beta1^-1) / (2 * N *T) * (  norminv(Omega)/(1-2*tau) )^2;
    epsilon = 2 * log(2/gamma -1);
    [Tardos_code,p]    = tardos_fp_code(epsilon, beta1, c0,c1);
    CONSIS = [];
    for t = 1:3
        K = ceil( topK(t) * 156 );
        consistency = 0;
        for j = 1:c1
            LDP_R = LDP_dbSNP(R,epsilon);
            MR = FP_dbSNP_one_seed(LDP_R,secretKey,gamma,Tardos_code(j,:));
            p_values = snp_p_value(MR,control_index, case_index);
            [~,idx] = sort(p_values,'descend');
            consistency = consistency + length( intersect (idx([1:K]),  ...
                ground_truth_idx([1:K])     )      )/K;
        end
        CONSIS = [CONSIS; consistency/c1];
    end
    CONSIS_two_stage = [CONSIS_two_stage CONSIS];
end

