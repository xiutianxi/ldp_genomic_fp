




clear;clc;close all;
addpath('saved data')
addpath('funcs')

load 'All_Pop.mat'

R = All_Pop;

secretKey = 'dbSNP';

sp_id = 10;


% random bit flipping versus per_cmp

fp = sp_id_fingerprint_generate(secretKey, sp_id);



epsilon = 5;

% MR = LDP_FP_dbSNP_one_seed(R,secretKey,epsilon,fp);

% MR = DPFP(R,secretKey,epsilon1,fp); switch to one seed to speed up, since
% only flipping attack is considered, use DPFP(...) if subset, superset
% attacks are considered, todo: efficient implementation of random sequnce
% under multiple seeds.. 

per_list = [10:5:45] / 100;

per_cmp = [];

for j = 1:length(per_list)
    per  = per_list(j);

    R_attack = random_bit_flipping(MR,per);

    f_extract = LDP_FP_dbSNP_extraction_one_seed(R_attack,R,secretKey,epsilon);

    f_extract = double(f_extract);
    per_cmp = [ per_cmp length(find(f_extract~=fp))/128 * 100];
end

per_cmp


