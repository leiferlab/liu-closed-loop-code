function no_weight_tuples = get_behavior_tuples(number_of_behaviors, degree)
%this generates all possible transition tuples

    iteration = 1;
    no_weight_tuples = 1:number_of_behaviors;
    while iteration < degree
        no_weight_tuples = combvec(no_weight_tuples, 1:number_of_behaviors);
        iteration = iteration + 1;
    end
    no_weight_tuples = no_weight_tuples';
end

