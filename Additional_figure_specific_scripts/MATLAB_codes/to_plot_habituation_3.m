clear
clc
% close all

load('/projects/LEIFER/Sandeep/APIData/20210624_RunFullWormRails_Sandeep_AML67_10ulRet_red/habituation_num_reversal_rails_AML67_20210624_3.mat')
load('/projects/LEIFER/Sandeep/APIData/20210614_RunFullWormRails_Sandeep_AML67_10ulRet_red/habituation_num_reversal_rails_AML67_20210614_3.mat')
load('/projects/LEIFER/Sandeep/APIData/20200902_RunFullWormRails_Sandeep_AML67_10ulret/habituation_num_reversal_rails_AML67_20200902_3.mat')

load('/projects/LEIFER/Sandeep/APIData/20210723_RunFullWormRails_Sandeep_AKS_483.7.e_mec4_Chrimson_10ulRet_red/habituation_num_reversal_rails_AML470_20210723_3.mat')
load('/projects/LEIFER/Sandeep/APIData/20210721_RunFullWormRails_Sandeep_AKS_483.7.e_mec4_Chrimson_10ulRet_red/habituation_num_reversal_rails_AML470_20210721_3.mat')
load('/projects/LEIFER/Sandeep/APIData/20210720_RunFullWormRails_Sandeep_AKS_483.7.e_mec4_Chrimson_10ulRet_red/habituation_num_reversal_rails_AML470_20210720_3.mat')

bins_of_interest=[0:3600:54000];

[mean_AML67,errhigh_AML67,errlow_AML67]=habituation_plot([habituation_num_reversal_rails_AML67_20200902_3(:,4:5);habituation_num_reversal_rails_AML67_20210614_3(:,4:5);habituation_num_reversal_rails_AML67_20210624_3(:,4:5)]...
    ,bins_of_interest);

[mean_AML470,errhigh_AML470,errlow_AML470]=habituation_plot([habituation_num_reversal_rails_AML470_20210720_3(:,4:5);habituation_num_reversal_rails_AML470_20210721_3(:,4:5);habituation_num_reversal_rails_AML470_20210723_3(:,4:5)]...
    ,bins_of_interest);
%%
x_data = [2:2:30];
figure;
plot(x_data,mean_AML470,'-sb','LineWidth',1.5)
hold on
err_AML470=errorbar(x_data,mean_AML470,errlow_AML470,errhigh_AML470,'color','b','LineWidth',1.5,'LineStyle', '-');
hold on;
plot(x_data,mean_AML67,'-or','LineWidth',1.5)
hold on
err_AML67=errorbar(x_data,mean_AML67,errlow_AML67,errhigh_AML67,'color','r','LineWidth',1.5); 
hold on;
xlabel('Time (s)') 
ylabel('Probability of reversal')
set(get(get(err_AML67,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(err_AML470,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('AML470','AML67')
ylim([0 0.9])
xlim([0 32])
yticks([0:0.3:0.9])
ax = gca;
legend box off
ax.FontSize = 14;
box off;
grid off;
%%
function [mean_of_binary_data_array,errhigh_array,errlow_array]=habituation_plot(data,bins_of_interest)
    
    defined_bins=discretize(data(:,2),bins_of_interest);
    final_data=[data defined_bins];

    mean_of_binary_data_array=[];
    errhigh_array=[];
    errlow_array=[];
    
    for mn=1:size(bins_of_interest,2)-1
        binned_data=final_data(final_data(:,3)==mn,:);
        [mean,errhigh,errlow] = bootstrap_mean_and_ci(10000,0.05,binned_data(:,1));
        
        mean_of_binary_data_array(mn,1)=mean;
        errhigh_array(mn,1)=errhigh;
        errlow_array(mn,1)=errlow;
        
    end
    
end


