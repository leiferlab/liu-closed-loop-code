function [mean_of_binary_data,errhigh,errlow] = bootstrap_mean_and_ci(nboot,alpha,binary_data)

    %%%% this function will enable us to find the mean and error bars (95% CI using bootstrap)

    ci_cutoff = round([nboot*alpha/2 nboot*(1-alpha/2)]');
    
    bootstrap_data = bootstrp(nboot,@mean,binary_data);   
    sorted_bootstrap_data = sort(bootstrap_data);
    ci_data = [sorted_bootstrap_data(ci_cutoff(1)), sorted_bootstrap_data(ci_cutoff(2))];

    mean_of_binary_data=sum(binary_data)/size(binary_data,1);
    errhigh = ((ci_data(2)-mean_of_binary_data));
    errlow  = ((mean_of_binary_data-ci_data(1)));

% % % % % %         ci_cutoff = round([nboot*alpha/2 nboot*(1-alpha/2)]');
% % % % % %         TRN_forward_stim_ON_bootstrap = bootstrp(nboot,@mean,binary_data);   %%%% 100 bootstraps for stim ON
% % % % % %         
% % % % % %         mean_TRN_reversal_rate_while_forward_ON=sum(binary_data)/size(binary_data,1);
% % % % % %         
% % % % % %         TRN_forward_reversal_multiples_of_baseline_ON = TRN_forward_stim_ON_bootstrap ;   %%%%% dividing the above two arrays to find the ratio
% % % % % %         mean_TRN_fullrails_reversal_multiples_of_baseline_ON=mean_TRN_reversal_rate_while_forward_ON; %%% to find the mean. Index 2=stim ON, index 1=stim OFF
% % % % % % 
% % % % % %         sorted_ratio_matrix_TRN_forward_ON = sort(TRN_forward_reversal_multiples_of_baseline_ON);
% % % % % %         ci_TRN_forward_ON = [sorted_ratio_matrix_TRN_forward_ON(ci_cutoff(1)), sorted_ratio_matrix_TRN_forward_ON(ci_cutoff(2))];
% % % % % %         
% % % % % %         mean_of_binary_data = [mean_TRN_fullrails_reversal_multiples_of_baseline_ON];
% % % % % % 
% % % % % %         errhigh = [(ci_TRN_forward_ON(2)-mean_TRN_fullrails_reversal_multiples_of_baseline_ON)];
% % % % % % 
% % % % % %         errlow  = [(mean_TRN_fullrails_reversal_multiples_of_baseline_ON-ci_TRN_forward_ON(1))];


end