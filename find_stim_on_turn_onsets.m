function true_valid_stim_final=find_stim_on_turn_onsets(stims_delivered_during_mochis_turns,er_turn_onset_indices)
%%% this function make sure that we are only selecting stim which occurs after the worm crossed ellipse ratio threshold   
    true_valid_stim=[];
    for xx=1:length(stims_delivered_during_mochis_turns)
        for yy=1:length(er_turn_onset_indices)
           dummy_pair= [stims_delivered_during_mochis_turns(xx,1) er_turn_onset_indices(yy,1)];
           if (dummy_pair(1)-dummy_pair(2))<=10 && (dummy_pair(1)-dummy_pair(2))>0
               true_valid_stim= [true_valid_stim dummy_pair(1)];
           end
        end
    end

    true_valid_stim_final=unique(true_valid_stim);
end
