function phi_dt = worm_phase_velocity(ProjectedEigenValues, parameters)
%This function outputs the phase velocity time series given the PCA time series% not used in paper

    phase_determining_PCs = [2, 3]; %which PCs are the sine and cosine of worm locomotion
    phase_determining_eigen_values = ProjectedEigenValues(phase_determining_PCs, :);
    
    x = phase_determining_eigen_values(1,:);
    y = phase_determining_eigen_values(2,:);

    %normalize x and y to be on the same scale (this is done in the paper)
    %the RMS of x and y are both 1
    x = x ./ parameters.PCxScale;
    y = y ./ parameters.PCyScale;
    y(y == 0) = eps;     % Avoid division by zero in phase calculation
    
    %define phase as invtan(x/y)
    phi = atan(x./y);
    NegYIndexes = find(y < 0);
    Index1 = find(phi(NegYIndexes) <= 0);
    Index2 = find(phi(NegYIndexes) > 0);
    phi(NegYIndexes(Index1)) = phi(NegYIndexes(Index1)) + pi;
    phi(NegYIndexes(Index2)) = phi(NegYIndexes(Index2)) - pi;
    
    %take the derivative off phi
    phi_dt = angdiff(phi(2:end),phi(1:end-1));
    
    %correct for discontinuity due to wrapping around and loss of a timepoint
    phi_dt = -[0, unwrap(phi_dt,[],1)]; 

    %gaussian smooth the result 
    phi_dt = smoothts(phi_dt, 'g', parameters.TrackingSmoothingWindow*parameters.SampleRate, parameters.TrackingSmoothingWindow*parameters.SampleRate);
    
%     %cap the phi_dt at some min and max
%     phi_dt(phi_dt < parameters.MinPhaseVelocity) = parameters.MinPhaseVelocity;
%     phi_dt(phi_dt > parameters.MaxPhaseVelocity) = parameters.MaxPhaseVelocity;
    
end

