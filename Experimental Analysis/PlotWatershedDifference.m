function [density_diff] = PlotWatershedDifference(embeddingValues1,embeddingValues2)
%Plots the density map along with the watershed; not used in paper
    maxVal = max(max(abs([embeddingValues1;embeddingValues2])));
    maxVal = round(maxVal * 1.1);
    load('reference_embedding.mat')

    sigma = 4; %change smoothing factor if necessary

    [xx,density1] = findPointDensity(embeddingValues1,sigma,501,[-maxVal maxVal]);
    [~,density2] = findPointDensity(embeddingValues2,sigma,501,[-maxVal maxVal]);
    
    density_diff = density1-density2;
    
    %get the average density difference
    L_flat = L(:);
    watershed_avg_density_diff = zeros(1,max(L_flat-1));
    density_diff_linear = density_diff(:); %linearize
    for watershed_region = 1:max(L(:)-1)
        watershed_avg_density_diff(watershed_region) = mean(density_diff_linear(L_flat == watershed_region));
    end

    max_avg_difference = max(abs(watershed_avg_density_diff));
    labeled_avg_diff = zeros(size(L_flat));
    
    for watershed_region = 1:max(L(:)-1)
        labeled_avg_diff(L_flat == watershed_region) = watershed_avg_density_diff(watershed_region);
    end
    labeled_avg_diff = reshape(labeled_avg_diff,size(L,1),size(L,2));

    maxDensity = max(abs(density_diff_linear(:)));
    [ii,jj] = find(encapsulate_watershed_matrix(L)==0);

    watershed_centroids = regionprops(L, 'centroid');
    watershed_centroids = vertcat(watershed_centroids.Centroid);
    watershed_centroids = round(watershed_centroids);

    %modify jet map
    my_colormap = redblue;

    %figure
    hold on
% %     %plot the point by point difference
% %     imagesc(xx,xx,density_diff)
% %     plot(xx(jj),xx(ii),'k.')
% %     axis equal tight off xy
% %     caxis([-maxDensity maxDensity])
    %plot the average difference
    imagesc(xx,xx,labeled_avg_diff)
    plot(xx(jj),xx(ii),'k.')
    axis equal tight off xy
    caxis([-max_avg_difference max_avg_difference])
    colormap(my_colormap)
%     for region_index = 1:size(watershed_centroids,1)-1
%         text(xx(watershed_centroids(region_index,1)), ...
%             xx(watershed_centroids(region_index,2)), ...
%             num2str(region_index), 'color', 'k', ...
%             'fontsize', 12, 'horizontalalignment', 'center', ...
%             'verticalalignment', 'middle');
%     end
    hold off
    colorbar
end

