function [result]= myNMIACCwithmean(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,size(U,2));
maxIter = 1;

indx = kmeans(U_normalized,numclass,'MaxIter',100, 'Replicates',3); %3
indx = indx(:);
result = Clustering8Measure(Y,indx); % result = [ACC nmi Purity Fscore Precision Recall AR Entropy];
ar =RandIndex(Y,indx);
end