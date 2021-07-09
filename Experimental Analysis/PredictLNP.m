function [predicted_behavior_rate] = PredictLNP(LEDPower, linear_kernel, modeled_fit, bin_size)
%Predicts the behavioral rate based on LNP model parameters; not used in paper

    if nargin < 4
        %no bin size specified
        bin_size = 1;
    end
    
    filtered_signal = padded_conv(LEDPower, linear_kernel);
    predicted_behavior_rate = modeled_fit(filtered_signal);
    predicted_behavior_rate = reshape(predicted_behavior_rate, [], bin_size);
    predicted_behavior_rate = mean(predicted_behavior_rate, 2)';
end

