

clear;clc;close all;
addpath('saved data')
addpath('funcs')

load 'All_Pop.mat'

R = All_Pop;
secretKey = 'dbSNP'; 



c0 = 3;

c1 = c0;


epsilon = ceil( 2*log(c0-1) );

beta1 = 10^-3;

[Tardos_code,p] = tardos_fp_code(epsilon, beta1, c0,c1);

MR_list = {};

for i = 1:c1
    MR_list{i} = LDP_FP_dbSNP_one_seed(R,secretKey,epsilon,Tardos_code(i,:));
end


% MR_content_list = {};

for i = 1:c1
    MR_content_list(:,:,i) =  MR_list{i}(:,[2:end]).Variables;
end

MR_content_majority = mode(MR_content_list,3);

R_collusion_attack = R;
R_collusion_attack(:,[2:end]).Variables = MR_content_majority;

L = size(Tardos_code,2);

f_extract = LDP_FP_dbSNP_extraction_one_seed(R_collusion_attack,R,secretKey,epsilon,L);

f_extract = double(f_extract);


Z  = 20 * c0 * (log(1/beta1));

accusation_sum = zeros(1,c1);

for i = 1:c1
    u = Tardos_code(i,:)*sqrt(  (1-p)/p ) + (  Tardos_code(i,:) -1 ) * sqrt(  p/(1-p) );
    accusation_sum(i) =  dot(f_extract, u);
end

sum(accusation_sum>Z)