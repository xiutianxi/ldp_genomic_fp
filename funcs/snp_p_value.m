function p_values = snp_p_value(R,control_index, case_index)


content = R(:,[2:end]).Variables;


control_content = content(control_index,:);

case_content = content(case_index,:);


snp_number = size(content,2);

p_values = [];

for i = 1:snp_number
    c0 = sum(control_content(:,i)==0);
    c1 = sum(control_content(:,i)==1);
    c2 = sum(control_content(:,i)==2);

    s0 = sum(case_content(:,i)==0);
    s1 = sum(case_content(:,i)==1);
    s2 = sum(case_content(:,i)==2);

    OR = c0 * (s1+s2) / (s0 * (c1+c2));

    stdErr = sqrt(   1/(s1+s2) + 1/s0 + 1/(c1+c2) + 1/c0     );

    z = log(OR) / stdErr;

    p = 1- normcdf(z) + normcdf(-z);

    p_values = [p_values p];

end