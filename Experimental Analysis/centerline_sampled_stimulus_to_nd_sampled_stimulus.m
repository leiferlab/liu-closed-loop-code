function nd_sampled_stimulus = centerline_sampled_stimulus_to_nd_sampled_stimulus(centerline_sampled_stimulus, nd)
    %nd is the number of "zones" to resample to
    n_centerline_pts = size(centerline_sampled_stimulus,2);
    spacings = linspace(1, n_centerline_pts+1, nd+1);
    nd_sampled_stimulus = zeros(size(centerline_sampled_stimulus,1), nd);
    for index = 1:nd
        start_pt = round(spacings(index));
        end_pt = round(spacings(index+1))-1;
        nd_sampled_stimulus(:,index) = no_mean_median(centerline_sampled_stimulus(:,start_pt:end_pt),2);
    end
end