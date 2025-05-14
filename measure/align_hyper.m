function [M_result]= align_hyper(M,RWofCluster1,alpha)
[a,~] = size(RWofCluster1);
M_result = zeros(size(M));
    for i = 1:a
        idx = RWofCluster1(:,i)>=alpha;
        temp = sum(M(:, idx),2);
        temp(temp>0) = 1;
        if sum(temp) == 0
        temp = M(:,i);
        end
        M_result(:,i) = temp;
        clear idx
    end
end

