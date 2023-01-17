function [f,f_as_0,f_as_1] = LDP_FP_dbSNP_extraction_one_seed(MR,R,secretKey,epsilon,L)



f_as_0 = zeros(L,1);

f_as_1 = zeros(L,1);


K = 2; % snp data has 2 bits

p = 1/( exp(epsilon/2)+1 );

[row_num,col_num] = size(R);


%% original bits

content  = ( R(:,[2:end]).Variables )';  % col_num by row_num

content_vec = content(:);


content_vec_binary = dec2bin(content_vec)-'0';


%% marked bits

content_marked  = ( MR(:,[2:end]).Variables )';  % col_num by row_num

content_marked_vec = content_marked(:);


content_marked_vec_binary = dec2bin(content_marked_vec)-'0';


%% recover B

B = xor(content_vec_binary, content_marked_vec_binary);

%%

rng(sum(double(double(secretKey))))

mark_bit_pos = rand(size(content_vec_binary))<p;
temp = content_vec_binary;

Dic_key_range = 2^25; % attribute number  = col_num-1

for k = 1:K

    rnd_seq = datasample([1:Dic_key_range],length(content_vec)*2,'Replace',false);

    rnd_seq = reshape(rnd_seq,length(content_vec),[]);




    B = xor( content_vec_binary(:,end-k+1) , content_marked_vec_binary(:,end-k+1) );


    x_bit = mod(  rnd_seq(:,1),2 );


    f_bit = xor(x_bit,B);

    l = mod(rnd_seq(:,2),L)+1;
    for j = 1:L
        f_as_0(j) = f_as_0(j) + sum(f_bit(l==j)==0);
        f_as_1(j) = f_as_1(j) + sum(f_bit(l==j)==1);
    end
end


f = double( (f_as_1>f_as_0)' );


end
