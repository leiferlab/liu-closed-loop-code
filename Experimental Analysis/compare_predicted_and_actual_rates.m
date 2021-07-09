function [ score ] = compare_predicted_and_actual_rates(predicted_rates, observed_events)
%computes a metric to see how predictive our model is against the actual
%observed rates; not used in paper
    fps = 14;
    dt = 1/(fps*60);
    %observed_rates = double(observed_rates)*fps*60; % put it in transitions/min
    rdt = predicted_rates.*dt;
    predicted_probability = rdt; %Bernouli process
    %predicted_probability = rdt.*exp(rdt); %poisson probability for one event (k=1)
    score = mean((predicted_probability.*observed_events)); %average squared residual
end

