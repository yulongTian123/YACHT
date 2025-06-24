%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
% This is a demo for the LWEA and LWGP algorithms. If you find this %
% code useful for your research, please cite the paper below.       %
%                                                                   %
% Dong Huang, Chang-Dong Wang, and Jian-Huang Lai.                  %
% "Locally weighted ensemble clustering."                           %
% IEEE Transactions on Cybernetics, 2018, 48(5), pp.1460-1473.      %
%                                                                   %
% The code has been tested in Matlab R2014a and Matlab R2015a on a  %
% workstation with Windows Server 2008 R2 64-bit.                   %
%                                                                   %
% https://www.researchgate.net/publication/316681928                %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ECI = computeECI_hyper(baseClsSegs, para_theta,n)

%M = size(mClsLabels,2);
ETs = getAllClsEntropy(baseClsSegs);
ECI = exp(-ETs./para_theta./n);
end

function Es = getAllClsEntropy(baseClsSegs)
% Get the entropy of each cluster w.r.t. the ensemble


[~, nCls] = size(baseClsSegs);

Es = zeros(nCls,1);
for i = 1:nCls
    temp = baseClsSegs(:,i);
    idx = temp~=0;
    set = baseClsSegs(idx,:);
    set(:, all(set==0)) = [];
    Es(i) = getOneClsEntropy(set);
    clear set
end
end
function E = getOneClsEntropy(set)
% Get the entropy of one cluster w.r.t the ensemble

% The total entropy of a cluster is computed as the sum of its entropy
% w.r.t. all base clusterings.
E = 0;
count = sum(set); 
if size(set,1) == 1
    return 
end

count = count./size(set,1);
E = E-sum(count.*log2(count));
end

