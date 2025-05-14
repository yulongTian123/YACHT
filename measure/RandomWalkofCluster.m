function R = RandomWalkofCluster(W)

    N = size(W,1);
    k = 10;
    beta = 1;

    W = W - diag(diag(W));
    D = diag(1 ./ sum(W, 1));
    D(isinf(D)) = 0;
    W_tilde = D * W;
    
    tmpO = W_tilde;
    O_tilde = W_tilde*W_tilde';
    
    for i = 1:k-1
        tmpO = tmpO*W_tilde;
        O_tilde = O_tilde + beta * (tmpO*tmpO');
    end
    O_i = repmat(diag(O_tilde), 1, N);
    R = O_tilde./sqrt(O_i.*O_i');
    
    sum_W = sum(W_tilde, 2);
    isolatedIdx = find(sum_W<10e-10);
    if numel(isolatedIdx)>0
        R(isolatedIdx,:) = 0;
        R(:,isolatedIdx) = 0;
    end
end