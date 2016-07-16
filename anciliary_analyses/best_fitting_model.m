function [  ] = best_fitting_model()

    import rsa.*
    import rsa.meg.*

    models_to_chose_from = { ...
        'FBK', ...
        'L2', ...
        'L3', ...
        'L6', ...
        'BN7', ...
        'triphone', ...
        'feature'};

    maps_base_path = '/imaging/cw04/CSLB/Lexpro/Analysis_DNN/CWD_win25_language_10242/';
    all_vals_template_template = fullfile(maps_base_path, 'Maps_%s/lexpro-bn-sl_group_t-map_tfce-%sh.stc');
    
    for chi = 'LR'
        for model_i = 1:numel(models_to_chose_from)
           model = models_to_chose_from{model_i};
           
           %% Load model and insert into stack
           
           this_model_path = sprintf(all_vals_template_template, model, lower(chi));
           
           stc_metadata = mne_read_stc_file1(this_model_path);
           
           if model_i == 1
               [n_vertices, n_timepoints] = size(stc_metadata.data);
               all_model_stack = zeros(n_vertices, n_timepoints, numel(models_to_chose_from));
           end
           
           all_model_stack(:, :, model_i) = stc_metadata.data;
           
        end
       
        %% Pick best models
        
        % Best models at each vertex,timepoint
        [max_vals, max_val_is] = max(all_model_stack, [], 3);
        
        % Best model over time-averaged values at each vertex
        [max_vals_overall, max_val_overall_is] = max(squeeze(mean(all_model_stack, 2)), [], 2);

        %% Write out maps
        
        max_vals_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'best_model_each_timepoint-%sh.stc'), lower(chi));
        max_vals_overall_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'best_model_overall-%sh.stc'), lower(chi));
        
        write_stc_file(stc_metadata, max_val_is, max_vals_path);
        write_stc_snapshot(stc_metadata, max_val_overall_is, max_vals_overall_path);
        
        %% Write out individual model-masked maps
        for model_i = 1:numel(models)
           masked_vals = zeros(size(max_val_is));
           masked_vals(max_val_is == model_i) = 1;
           
           masked_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'model_%d_%s_masked-%sh.stc'), model_i, model, lower(chi));
           write_stc_file(stc_metadata, masked_vals, masked_path);
        end
        
    end
    
end

