function MSE = compare_predicted_and_actual_rates_MSE(predicted_rates, observed_events)
%computes the MSE to see how predictive our model is against the actual
%observed rates using poisson distribution; not used in paper
    fps = 14;
    dt = 1/(fps*60);
    %observed_rates = double(observed_rates)*fps*60; % put it in transitions/min
    rdt = predicted_rates.*dt;
    predicted_probability = rdt; %Bernouli process
    %predicted_probability = rdt.*exp(rdt); %poisson probability for one event (k=1)
    MSE = mean((predicted_probability-observed_events).^2); %average squared residual
%    MSE = mean((predicted_probability.*observed_events)); %average squared residual
end