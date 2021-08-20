%%%% This file generates the figure 4b in the publication: "Closed-loop targeted optogenetic stimulation of C. elegans populations"

clear
clc
close all

%%%%%% loading closed loop data from 20210723
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210723_turns_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210723_turns_80uW_3s_new4.mat')

%%%%%% loading closed loop data from 20210721
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210721_turns_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210721_turns_80uW_3s_new4.mat')

%%%%%% loading closed loop data from 20210720
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210720_turns_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/closed_loop/AML470_20210720_turns_80uW_3s_new4.mat')

%%%%%% loading open loop data from 20210723
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210723_rails_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210723_rails_80uW_3s_new4.mat')

%%%%%% loading open loop data from 20210721
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210721_rails_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210721_rails_80uW_3s_new4.mat')

%%%%%% loading open loop data from 20210720
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210720_rails_0uW_3s_new4.mat')
load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/Matlab_codes_and_data/data_to_plot_probability_of_reversal/AML470/open_loop/AML470_20210720_rails_80uW_3s_new4.mat')

%% combining data from all the three days
%%%% 80 uW closed loop data
turning_data_80uW=[AML470_20210723_turns_80uW_3s_new4
    AML470_20210721_turns_80uW_3s_new4
    AML470_20210720_turns_80uW_3s_new4];

%%%% 0.5 uW closed data
turning_data_0uW=[AML470_20210723_turns_0uW_3s_new4
    AML470_20210721_turns_0uW_3s_new4
    AML470_20210720_turns_0uW_3s_new4];

%%%%%%%%%

%%%% 80 uW open loop data
rails_data_80uW=[AML470_20210723_rails_80uW_3s_new4
    AML470_20210721_rails_80uW_3s_new4
    AML470_20210720_rails_80uW_3s_new4];

%%%% 0.5 uW open loop data
rails_data_0uW=[AML470_20210723_rails_0uW_3s_new4
    AML470_20210721_rails_0uW_3s_new4
    AML470_20210720_rails_0uW_3s_new4];

%% finding the mean and 95% ci using bootstrap function (bootstrap_mean_and_ci.m)

[mean_data(:,1),err_high(:,1),err_low(:,1)]=bootstrap_mean_and_ci(10000,0.05,rails_data_80uW(:,4));
[mean_data(:,2),err_high(:,2),err_low(:,2)]=bootstrap_mean_and_ci(10000,0.05,turning_data_80uW(:,4));
[mean_data(:,3),err_high(:,3),err_low(:,3)]=bootstrap_mean_and_ci(10000,0.05,rails_data_0uW(:,4));
[mean_data(:,4),err_high(:,4),err_low(:,4)]=bootstrap_mean_and_ci(10000,0.05,turning_data_0uW(:,4));

x = [1:2 4:5];
figure1=figure;
b=bar(x,mean_data,'FaceColor','flat');                

hold on

er = errorbar(x,mean_data,err_low,err_high);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth=1.5;
b.CData(1,:) = [0.07 0.62 1];
b.CData(2,:) = [1 0.41 0.16];
b.CData(3,:) = [0.07 0.62 1];
b.CData(4,:) = [1 0.41 0.16];

% % x_label_names={'  Fwd \newline(Expt.)','  Turn \newline (Expt.)','    Fwd \newline  (Cont.)','    Turn \newline  (Cont.)'};
x_label_names={'Fwd','Turn','Fwd','Turn'};

set(gca,'xticklabel',x_label_names,'FontSize',16)
yticks([0:0.2:0.8])
ylim([0 0.7])
ylabel('Probability of reversal')
box off
hold off

% Create textbox
annotation(figure1,'textbox',...
    [0.460571428571429 0.7 0.182285714285714 0.0937619047619057],...
    'String','AML470',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.0212857142857142 0.859523809523812 0.0590714285714286 0.0890000000000021],...
    'String','b.',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create line
annotation(figure1,'line',[0.28 0.40],...
    [0.845 0.845],'LineWidth',1.5);

% Create textbox
annotation(figure1,'textbox',...
    [0.300 0.845 0.089 0.073],...
    'String','***',...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create line
annotation(figure1,'line',[0.64 0.76],[0.239 0.239],...
    'LineWidth',1.5);

% Create textbox
annotation(figure1,'textbox',...
    [0.660 0.239 0.089 0.073],...
    'String',{'*'},...
    'HorizontalAlignment','center',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create line
annotation(figure1,'line',[0.591071428571428 0.805357142857142],...
    [0.0323333333333337 0.0333333333333337]);

% Create line
annotation(figure1,'line',[0.230357142857143 0.444642857142857],...
    [0.0323333333333339 0.0333333333333339]);

% % % % % % Create textbox
% % % % % annotation(figure1,'textbox',...
% % % % %     [0.597499999999998 0.000357142857142849 0.204285714285716 0.073904761904763],...
% % % % %     'String','Control',...
% % % % %     'HorizontalAlignment','center',...
% % % % %     'FontSize',14,...
% % % % %     'FitBoxToText','off',...
% % % % %     'EdgeColor','none');
% % % % % 
% % % % % % Create textbox
% % % % % annotation(figure1,'textbox',...
% % % % %     [0.233214285714284 0.000357142857142858 0.204285714285716 0.0739047619047629],...
% % % % %     'String','Experiment',...
% % % % %     'HorizontalAlignment','center',...
% % % % %     'FontSize',14,...
% % % % %     'FitBoxToText','off',...
% % % % %     'EdgeColor','none');

%%%%% code for doing statistical test
%%%%% for bar 1 and bar 2
n1_expt=size(rails_data_80uW,1);
n2_expt=size(turning_data_80uW,1);
x1_bar_expt=mean_data(:,1);
x2_bar_expt=mean_data(:,2);

p_hat_expt=(n1_expt*x1_bar_expt+n2_expt*x2_bar_expt)/(n1_expt+n2_expt);
t_value_expt=(x1_bar_expt -x2_bar_expt)/sqrt((p_hat_expt*(1-p_hat_expt)*(1/n1_expt+1/n2_expt)));
df_expt=n1_expt+n2_expt-2;
p_expt=2*(1-tcdf(abs(t_value_expt),df_expt))

%%%%% for bar 3 and bar 4
n1_control=size(rails_data_0uW,1);
n2_control=size(turning_data_0uW,1);
x1_bar_control=mean_data(:,3);
x2_bar_control=mean_data(:,4);

p_hat_control=(n1_control*x1_bar_control+n2_control*x2_bar_control)/(n1_control+n2_control);
t_value_control=(x1_bar_control -x2_bar_control)/sqrt((p_hat_control*(1-p_hat_control)*(1/n1_control+1/n2_control)));
df_control=n1_control+n2_control-2;
p_control=2*(1-tcdf(abs(t_value_control),df_control))
