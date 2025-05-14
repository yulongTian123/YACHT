function YY_label = base_clustering_preserve(mClsLabels,K)
[table] = tabulate(mClsLabels(:,2)) ;                 
[~,idy] = sort(table(:,2),'descend');
[a,~] = size(mClsLabels);
YY_label = zeros(a,K);
for i = 1:K
x = find(mClsLabels(:,2) ==idy(i));
YY_label(mClsLabels(x,1),i) = 1;
end


end