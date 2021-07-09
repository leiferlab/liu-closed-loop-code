function M = moving_maximum(A,k,use_absolute_max)
%calculates the moving max for A with sliding window k
    k_half = floor(k/2);
    M = zeros(size(A));
    for i = 1:size(A,2)
       if use_absolute_max
           [M(:,i), max_indecies] = max(abs(A(max(i-k_half,1):min(i+k_half,size(A,2)))),[],2);
           signs = sign(A(max(i-k_half,1):min(i+k_half,size(A,2))));
           M(:,i) = M(:,i) .* signs(:,max_indecies);
       else
           M(:,i) = max(A(max(i-k_half,1):min(i+k_half,size(A,2))),[],2);
       end
    end
end

