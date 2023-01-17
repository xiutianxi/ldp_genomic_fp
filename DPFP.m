function MR = DPFP(R,secretKey,epsilon1,fp)


K = 2;

p = 1/( exp(epsilon1/K)+1 );


bits_required =  floor( log2( max( R.Variables,[],1) ) )+1;
bits_required = bits_required(2:end);




[row_num,col_num] = size(R);

L = length(fp);
Start = 1;
% Stop = row_num;

T = col_num-1;

Dic_key_range = row_num*K*T; % attribute number  = col_num-1




for i  = 1:row_num
    
    if mod(i,100)==0
        i
    end
% i
    primary_key_att_value = R{i,1};
    for t = 1:T
        bit_length = bits_required(t);
        
        k_min = min(bit_length,K);
        
        r_it = R{i,t+1}; % skip the first column of index, i.e., primary key
        r_it_binary =    dec2bin(r_it);
        temp = int2str( zeros( 1, k_min-length(r_it_binary)  ) );
        temp = temp(find(~isspace(temp)));
        r_it_binary = [  temp   r_it_binary]; % zero padding
        
        
        for k = 1:k_min
%                         [i t k]
            seed = [double(secretKey)  primary_key_att_value t k];
            rng(sum(seed));
            rnd_seq = datasample([1:Dic_key_range],3,'Replace',false);
            %             if ~mod(rnd_seq(1),floor( 1/(2*p) ) )
            if rand<2*p
                
                x = mod(  rnd_seq(2),2 );
                
                l = mod(rnd_seq(3),L)+1;
                
                f = fp(l);
                
                B = xor(x,f);
                
                
                r_it_binary(end-k+1) = int2str(  xor(B,str2num(r_it_binary(end-k+1)))  );
            end
        end
        r_it= bin2dec(r_it_binary);
        
        R{i,t+1} = r_it;
        %         i
    end
    
end

MR = R;

end
