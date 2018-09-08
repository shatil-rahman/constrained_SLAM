function edges = edge_counter(l,r)

    edge_counts = zeros(17,1);
    for i=1:length(l)
        edge_counts(i) = nnz(r(:,i));
    end

    [b,i] = sort(edge_counts);
    edges = nchoosek(i(end-3:end),2);
    edges = edges(1:end-1,:);
    
end
