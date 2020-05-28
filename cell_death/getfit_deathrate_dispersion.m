function sumsq_error = getfit_deathrate_dispersion(log_alpha_dispersion, moi_vector_3hpi, moi_vector, data_allhpi, FACs_data_18hpi, hpi, model_num, which_distribution)

alpha_dispersion = exp(log_alpha_dispersion);

% only fit to the first 4 treatment groups of FACs data
these_moi_FACs = moi_vector(1:(end-3));

% switch-case over models: 1-4
switch model_num
    
    case 1
        
        % 'model 1': constant rate of cell death, input-independent
        alpha = alpha_dispersion(1);
        k_dispersion = alpha_dispersion(end);
        
    case 2
        % 'model 2': Weibull hazard function cell death rate, input-independent
        alpha = alpha_dispersion(1:2);
        k_dispersion = alpha_dispersion(end);
        
    case 3
        % 'model 3': constant rate of cell death (time-independent), input-dependent
        alpha = alpha_dispersion(1:2);
        k_dispersion = alpha_dispersion(end);
        
    case 4
        % 'model 4': time-dependent + input-dependent
        alpha = alpha_dispersion(1:3);
        k_dispersion = alpha_dispersion(end);
        
end



sumsq_error = 0;
for n = 1:length(hpi)
    
    % this hour post infection
    this_hpi = hpi(n);
    
    % data at this hpi
    data_this_hpi = data_allhpi(n,:);
    
    if n == 1    
        these_moi = moi_vector_3hpi;
    else
        these_moi = moi_vector;
    end
    
    
    if n == 4 % 18 hpi
        
        percent_alive = get_percentalive(alpha,this_hpi,these_moi,k_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((data_this_hpi - percent_alive).^2);
        
        %         sumsq_error_percentalive_hpi(counter) = sum((data_this_hpi - percent_alive).^2);
        
        percent_aliveinfected = get_percentaliveinfected(alpha,this_hpi,these_moi_FACs,k_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((FACs_data_18hpi - percent_aliveinfected).^2);
        
        %         sumsq_error_total_m2 = sumsq_error;
        %         sumsq_error_percentaliveinfected = sumsq_error_total_m2 - sum(sumsq_error_percentalive_hpi_m2);
        
    else
        percent_alive = get_percentalive(alpha,this_hpi,these_moi,k_dispersion,model_num,which_distribution);
        sumsq_error = sumsq_error + sum((data_this_hpi - percent_alive).^2);
        
        %         sumsq_error_percentalive_hpi(counter) = sum((data_this_hpi - percent_alive).^2);
        %         counter = counter + 1;
        
    end
    
    
    
end

