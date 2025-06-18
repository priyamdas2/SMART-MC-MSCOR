function T = extract_dynamic_effect(M_dummy, N, var_names, row_labels_initials, Trt_trans, row_labels)

    % Extract initial state probabilities
    output_matrix_initials = zeros(N, 100);
    for i = 1:N
        output_matrix_initials(i, :) = squeeze(M_dummy(1, i, :));
    end
    T_initials = array2table(output_matrix_initials, 'VariableNames', var_names(2:end));
    T_initials = addvars(T_initials, row_labels_initials', 'Before', 1, 'NewVariableNames', 'Transition');

    % Extract transition probabilities for selected (u,v) pairs
    output_matrix_trans = zeros(size(Trt_trans,1), 100);
    for i = 1:size(Trt_trans,1)
        u = Trt_trans(i, 1) + 1; % convert u from 0-based to 1-based index
        v = Trt_trans(i, 2);
        output_matrix_trans(i, :) = squeeze(M_dummy(u, v, :));
    end
    T_trans = array2table(output_matrix_trans, 'VariableNames', var_names(2:end));
    T_trans = addvars(T_trans, row_labels', 'Before', 1, 'NewVariableNames', 'Transition');

    % Combine both
    T = [T_initials; T_trans];

end
