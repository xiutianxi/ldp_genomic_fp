function LDP_R = LDP_dbSNP(R,epsilon)

content = R(:,[2:end]).Variables;

p = 1/( exp(epsilon) + 2 );


index = rand(size(content))<p;

content_tobe_perturbed = content( index );


% define mapping 0-[1 2], 1-[0 2], 2-[0 1]


MAP = zeros(length(content_tobe_perturbed),2);

MAP(find(content_tobe_perturbed==0),:) = repmat([1 2],sum(content_tobe_perturbed==0),1);
MAP(find(content_tobe_perturbed==1),:) = repmat([0 2],sum(content_tobe_perturbed==1),1);
MAP(find(content_tobe_perturbed==2),:) = repmat([0 1],sum(content_tobe_perturbed==2),1);


noise_index = datasample([1,2],length(content_tobe_perturbed),'Weights',[0.5 0.5]);

noise_index = (noise_index-1) * length(content_tobe_perturbed) ...
    + [1:length(content_tobe_perturbed)];

ldp_noise = MAP(noise_index);

content(index) = ldp_noise;

LDP_R = R;
LDP_R(:,[2:end]).Variables = content;

end