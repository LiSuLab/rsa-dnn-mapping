function [  ] = best_fitting_model()

    import rsa.*
    import rsa.meg.*

    % TODO: put correct model names in here
    models_to_chose_from = { ...
        'FBK', ...
        'L2', ...
        'L3', ...
        'L6', ...
        'BN7', ...
        'triphone', ...
        'feature'};

    maps_base_path = '/imaging/cw04/CSLB/Lexpro/Analysis_NN_mapping/CWD_win25_language_10242/';
    all_vals_template_template = fullfile(maps_base_path, 'Maps_%s/lexpro-bn-sl_group_t-map_tfce-%sh.stc');
    
    for chi = 'LR'
        for model_i = 1:nume(models_to_chose_from)
           model = models_to_chose_from{model_i};
           
           %% Load model and insert into stack
           
           this_model_path = sprintf(all_vals_template_template, model, lower(chi));
           
           stc_metadata = mne_read_stc_file_1(this_model_path);
           
           if model == 1
               [n_vertices, n_timepoints] = size(stc_metadata.data);
               all_model_stack = zeros(n_vertices, n_timepoints, numel(models_to_chose_from));
           end
           
           all_model_stack(:, :, model_i) = stc_metadata;
           
        end
       
        %% Pick best models
        
        % Best models at each vertex,timepoint
        [max_vals, max_val_is] = max(all_model_stack, [], 3);
        
        % Best model over all time at each vertex
        [max_vals_overall, max_val_overall_ts] = max(max_vals, [], 2);
        max_val_overall_is = max_val_is(:, max_val_overall_ts);

        %% Write out maps
        
        max_vals_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'best_model_each_timepoint_%sh.stc'), lower(chi));
        max_vals_overall_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'best_model_overall_%sh.stc'), lower(chi));
        
        write_stc_file(stc_metadata, max_val_is, max_vals_path);
        write_stc_snapshot(stc_metadata, max_val_overall_is, max_vals_overall_path);
        
    end
    
end

