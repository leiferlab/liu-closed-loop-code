function discritized_v = discretize_given_edges(v, edges)
%similar to the matlab discritize function in later versions
    discritized_v = zeros(size(v));
    discritized_v(v < edges(1)) = 1;
    for edge_index = 1:numel(edges)-1
        discritized_v(v >= edges(edge_index) & v < edges(edge_index+1)) = edge_index;
    end
    discritized_v(v >= edges(end)) = numel(edges)-1;
end

