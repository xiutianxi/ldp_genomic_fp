function R_attack = random_bit_flipping(R,per)





[row_num,col_num] = size(R);

content  = ( R(:,[2:end]).Variables )';  % col_num by row_num

content_vec = content(:);


content_vec_binary = dec2bin(content_vec)-'0';



mark_bit_pos = rand(size(content_vec_binary))<per;

temp = content_vec_binary;

temp(mark_bit_pos==1) = 1-content_vec_binary(mark_bit_pos==1);


temp_bin2dec = bin2dec( char(temp+'0') );

content_mark =    reshape(  temp_bin2dec,  col_num-1, row_num)  ;


content_mark(content_mark>2)=0;

R_attack = R;
R_attack(:,[2:end]).Variables = content_mark';








end