function max_v = absolutevalue_maximum(v)
%calculates the maximums of the absolute value for A

    [max_v, max_index] = max(abs(v),[],2);
    signs = sign(v);
    max_v = max_v .* signs(max_index);
end

