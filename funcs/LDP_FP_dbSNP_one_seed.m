function MR = LDP_FP_dbSNP_one_seed(R,secretKey,epsilon,fp)

L = length(fp);

K = 2; % snp data has 2 bits

p = 1/( exp(epsilon/2)+1 );

[row_num,col_num] = size(R);

content  = ( R(:,[2:end]).Variables )';  % col_num by row_num

content_vec = content(:);


content_vec_binary = dec2bin(content_vec)-'0';

rng(sum(double(double(secretKey))))

mark_bit_pos = rand(size(content_vec_binary))<p;
temp = content_vec_binary;

Dic_key_range = 2^25; % attribute number  = col_num-1

for k = 1:K

    rnd_seq = datasample([1:Dic_key_range],length(content_vec)*2,'Replace',false);

    rnd_seq = reshape(rnd_seq,length(content_vec),[]);




    x_bit = mod(  rnd_seq(:,1),2 );

    l = mod(rnd_seq(:,2),L)+1;

    f_bit = fp(l);

    B = xor(x_bit,f_bit');

    temp(:,end-k+1) = xor(B, temp(:,end-k+1) );
end




temp(mark_bit_pos==0) = content_vec_binary(mark_bit_pos==0);


temp_bin2dec = bin2dec( char(temp+'0') );

content_mark =    reshape(  temp_bin2dec,  col_num-1, row_num)  ;


content_mark(content_mark>2)=0;

MR = R;
MR(:,[2:end]).Variables = content_mark';


end
