function [] = plot_worm_frame(I, center_lines, centerline_color, UncertainTips, centerline_stimulus, debugimage, additional_text)
%     IWFig = findobj('Tag', ['IWFig', num2str(plotting_index)]);
%     if isempty(IWFig)
%         IWFig = figure('Tag', ['IWFig', num2str(plotting_index)]);
%     else
%         figure(IWFig);
%     end
    %used for debugging
%     hold off;

%     clf
    if nargin > 5
        %debugimage inputted
        I = I * 2;
        imshow(cat(3, I, I, I) + debugimage, [], 'InitialMagnification', 500, 'Border','tight');
    else
        imshow(I, [], 'InitialMagnification', 500, 'Border','tight');
    end
    
    if isempty(centerline_color)
        centerline_color = [1, 1, 0];
    end
    

    hold on
    %head
    plot(center_lines(1,2), center_lines(1,1), '.', 'Color', centerline_color, 'markersize',30)
    plot(center_lines(:,2), center_lines(:,1), '-', 'Color', centerline_color, 'LineWidth',4)
%     if ~isempty(centerline_stimulus)
%         for centerline_point_index = 1:size(center_lines,1)
%             if isnan(centerline_stimulus(centerline_point_index))
%                 dot_color = [0 1 0];
%             elseif centerline_stimulus(centerline_point_index) > 0
%                 dot_color = [centerline_stimulus(centerline_point_index)/255, 0, 0];
%             elseif centerline_stimulus(centerline_point_index) < 0
%                 dot_color = [0,0,-centerline_stimulus(centerline_point_index)/255];
%             else
%                 dot_color = [0 0 0];
%             end
%             plot(center_lines(centerline_point_index,2), center_lines(centerline_point_index,1), '.', 'Color', dot_color, 'markersize',20)
%         end
%     end
    %plot(center_lines(1,2), center_lines(1,1), '.g','markersize',50)

    %uncertain tips
    if ~isempty(UncertainTips) && ~isempty(UncertainTips.Tips)
        plot(UncertainTips.Tips(:,2), UncertainTips.Tips(:,1), 'oy')
    end
%     %tail
%     plot(center_lines(end,2), center_lines(end,1), 'ob')
% 
%     
    %direction
%     quiver(size(I,2)/2, size(I,1)/2, sind(direction)*speed*100, -cosd(direction)*speed*100, 'AutoScale','off', 'Linewidth', 1.5);
%     %title (['Eccentricity = ', num2str(eccentricity)]);    
%     title (['( ', num2str(centroid(1,1)),' , ', num2str(centroid(1,1)),' )']);
    %score
    if ~isempty(additional_text)
        text(10, 10, additional_text, 'Color', 'y');
    end
%     %eccentricity
%     text(10, 60, num2str(eccentricity), 'Color', 'y');
    hold off;
%     drawnow
end