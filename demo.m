%You can adjust the \rho and \mu parameters in the solve_YACHT.m file 
% to achieve better results.
%%
clear
data = 'Segment_base_clustering.mat';
load(data);
%n = 20;
gt = y;
members = E;
n = 20;
k = length(unique(gt));
anchor_all = [k-1,k:k+3];
alpha_all = 0.9:0.01:1;
for ia = 1:length(anchor_all)   
    
for alpha_i = 1:length(alpha_all)
    alpha = alpha_all(alpha_i);
for i= 1:10
    [a,b] = size(members);
    zz = RandStream('mt19937ar','Seed',i);
    RandStream.setGlobalStream(zz);
    indx = randperm(b);
    EC_end = members(:,indx(1:n));% Base clustering result
    M = Gbe(EC_end);
    para_theta = 0.4;
    [~, mClsLabels] = computeMicroclusters(M);
    YY_label = base_clustering_preserve(mClsLabels,k);
    simOfCluster = full(simxjac(M'));
    RWofCluster = RandomWalkofCluster(simOfCluster);
    M_new = align_hyper(M,RWofCluster,alpha);
    consensus_ECI = computeECI_hyper(M_new, para_theta,n);
    [U,F] = solve_YACHT(M_new,anchor_all(ia),consensus_ECI',YY_label);
    res = myNMIACCwithmean(U,gt,k);
    ACC(i) = res(1);
    ARI(i) = res(7);
end

mean_acc(ia,alpha_i)= mean(ACC);
mean_ari(ia,alpha_i) = mean(ARI);
std_acc(ia,alpha_i)= std(ACC);
std_ari(ia,alpha_i) = std(ARI);
end
end
