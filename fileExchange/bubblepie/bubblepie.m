function bubblepie(xlist,ylist,slist,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
% A Bubble Plot, Pie Chart Combination
% bubblepie(xlist,ylist,slist,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
% 
% Creates a plot with pie charts at (xlist, ylist) using graph_data, having
% size of the pie scaled by slist.  Graph_labels contains a title for each
% pie chart (it can be empty [] to avoiid labeling), Graph_legend indicates
% the contents of graph_data, and lab is a binary value indicating whether
% pie chart labels are displayed. Pie colormaps can be modified afterwards
% with the colormap function or by acessing indivdual axes children.
% 
% Example:
% x = -pi:1:pi;
% x = x';
% y = sin(x);
% s = 1.1+cos(x);
% graph_labels = mat2cell(x,ones(1,length(x)),1);
% graph_data = 10*rand(length(x),3);
% graph_legend = {'one','two','three'};
% xlab = 'radians';
% ylab = 'sin(x)';
% lab = 1;
% 
% bubblepie(x,y,s,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
% 
% title('BubblePie Plot')
% 
% To avoid titling pie charts pass in empty vector:
% bubblepie(x,y,s,graph_data,[],graph_legend,xlab,ylab,lab)
% 
% To indicate no pie labels are wanted by setting the parameter lab = 0:
% bubblepie(x,y,s,graph_data,graph_labels,graph_legend,xlab,ylab,0)
% 
%   Abraham Anderson
%   July 30, 2007, updated Oct 27, 2017

%check for empty labels
if ~exist('graph_labels') || isempty(graph_labels)
    graph_labels = repmat({''},length(graph_data),1);
end

graph_max_size = 0.10;
graph_min_size = 0.05;
graph_range = graph_max_size-graph_min_size;

canvas_max = 1-graph_max_size/2;
canvas_min = 0.1;
canvas_range = canvas_max-canvas_min;

maxx = max(xlist);
maxy = max(ylist);
minx = min(xlist) - 1;
miny = min(ylist) - 1;
maxs = max(slist);

figure('position', [0 0 600 600])
h0 = axes('position',[canvas_min,canvas_min,canvas_range,canvas_range], ...
    'xlim',[minx maxx],'ylim',[miny maxy]);
set(get(gca,'XLabel'),'String',xlab)
set(get(gca,'YLabel'),'String',ylab)
% text(0.1*maxx+minx,maxy-0.1*maxy,{'PieChart Groups (CCW):' graph_legend{:}},'verticalalignment','top');
if lab
    for i = 1:size(graph_data,1)
        s = slist(i)/maxs*graph_range+graph_min_size;
        x = (xlist(i)-minx)/(maxx-minx)*canvas_range+canvas_min-s/2;
        y = (ylist(i)-miny)/(maxy-miny)*canvas_range+canvas_min-s/2;
        d = graph_data(i,:);
        if sum(d)==0
            continue
        end
        d(d==0) = 0.1;
        axes('position',[x y s s])
        pie(d)
        title(graph_labels{i})
    end
else
    for i = 1:size(graph_data,1)
        s = slist(i)/maxs*graph_range+graph_min_size;
        x = (xlist(i)-minx)/(maxx-minx)*canvas_range+canvas_min-s/2;
        y = (ylist(i)-miny)/(maxy-miny)*canvas_range+canvas_min-s/2;
        d = graph_data(i,:);
        if sum(d)==0
            continue
        end
        d(d==0) = 0.1;
        axes('position',[x y s s])
        pie(d,repmat({''},length(d),1))
        title(graph_labels{i})
    end
end

set(gcf,'Currentaxes',h0)
