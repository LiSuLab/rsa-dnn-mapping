function [  ] = best_fitting_model()

    import rsa.*
    import rsa.meg.*
    import rsa.util.*
    
    %                                p <  0.05   0.01   0.001   0.0001
    vertex_level_thresholds.FBK.L      = [386.6, 662.8, 1370.5, 1488.3];
    vertex_level_thresholds.FBK.R      = [402.0, 973.5, 1731.9, 1859.1];
    
    vertex_level_thresholds.L2.L       = [447.9, 746.9, 1502.3, 1526.0];
    vertex_level_thresholds.L2.R       = [489.1, 1107.6, 1654.8, 1691.2];
    
    vertex_level_thresholds.L3.L       = [390.8, 744.4, 1238.2, 1240.1];
    vertex_level_thresholds.L3.R       = [446.5, 1152.6, 1696.2, 1717.0];
    
    vertex_level_thresholds.L4.L       = [418.1, 723.6, 1161.6, 1197.4];
    vertex_level_thresholds.L4.R       = [472.0, 925.5, 1926.4, 2163.1];
    
    vertex_level_thresholds.L5.L       = [449.6, 721.9, 1132.9, 1185.3];
    vertex_level_thresholds.L5.R       = [505.0, 906.8, 2008.4, 2044.4];
    
    vertex_level_thresholds.L6.L       = [402.1, 662.6, 1137.6, 1251.2];
    vertex_level_thresholds.L6.R       = [491.3, 1058.8, 1852.4, 1951.3];
    
    vertex_level_thresholds.BN7.L      = [414.7, 732.7, 1068.3, 1110.5];
    vertex_level_thresholds.BN7.R      = [463.0, 852.3, 1860.2, 1860.2];
    
    %vertex_level_thresholds.triphone.L = [422.9, 727.3, 1435.6, 1719.4];
    %vertex_level_thresholds.triphone.R = [472.9, 685.1, 2041.1, 2253.7];
    
    % 1: p < 0.05
    % 2: p < 0.01
    % 3: p < 0.001
    % 4: p < 0.0001
    threshold_level = 3;
    

    models_to_chose_from = fieldnames(vertex_level_thresholds);

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
           
           data_mesh = stc_metadata.data;
           
           % Theshold
           vlt = vertex_level_thresholds.(model).(chi)(threshold_level);
           data_mesh(data_mesh < vlt) = 0;
           
           all_model_stack(:, :, model_i) = data_mesh;
           
        end
       
        %% Pick best models
        
        % There's probably a smart way to do this, but I'm just going to do
        % it in a dumb loop so I can get it done.
        
        max_val_is = zeros(n_vertices, n_timepoints);
        
        for v = 1:n_vertices
            for t = 1:n_timepoints
                model_fits = squeeze(all_model_stack(v, t, :));
                
                if sum(model_fits(:)) == 0
                    % No models fit here, so leave it zero
                else
                    [max_val, max_val_is(v, t)] = max(model_fits);
                end
            end
        end
        
        %% Write out maps
        
        max_vals_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'best_model-%sh.stc'), lower(chi));
        write_stc_file(stc_metadata, max_val_is, max_vals_path);
        
        %% Write out individual model-masked maps
        for model_i = 1:numel(models_to_chose_from)
            model = models_to_chose_from{model_i};
            
            model_masked_vals = zeros(size(max_val_is));
            model_masked_vals(max_val_is == model_i) = 1;
           
            masked_path = sprintf(fullfile(maps_base_path, 'Summary_maps', 'model_%s-%sh.stc'), model, lower(chi));
            write_stc_file(stc_metadata, model_masked_vals, masked_path);
            
            %% And display the numbers!
            
            % For all sig models
            extent_count = zeros(1, n_timepoints);
            for t = 1:n_timepoints
               vertices_this_timepoint = all_model_stack(:, t, model_i);
               extent_count(t) = sum(vertices_this_timepoint(:) > 0);
            end
            count_string = sprintf('%d, ', extent_count);
            prints('%s-h SIG model %d %s: [%s]', chi, model_i, model, count_string);
            
            % For the best models
            extent_count = zeros(1, n_timepoints);
            for t = 1:n_timepoints
                vertices_this_timepoint = model_masked_vals(:, t);
                extent_count(t) = sum(vertices_this_timepoint(:));
            end
            count_string = sprintf('%d, ', extent_count);
            prints('%s-h BEST model %d %s: [%s]', chi, model_i, model, count_string);
            
        end
    end
    
end

