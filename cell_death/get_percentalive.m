function percent_alive = get_percentalive(alpha, hpi, moi, k_dispersion, model_num, which_distribution)


i_vals = 0:1000;

switch model_num
    
    case 1
        
        % cell death rate: time-independent, input-independent - a
        a = alpha(1);
        int_psi = a*hpi; % integrating
        
    case 2
        
        % cell death rate: time-dependent, input-independent - b*k*hpi^(k-1)
        k = alpha(1); % shape parameter
        b = alpha(2); % scale parameters
        int_psi = b*hpi^k; % integrating Weibull hazard function
        
    case 3
        
        % cell death rate: time-independent, input-dependent - c*i^eps
        c = alpha(1);
        eps = alpha(2);
        int_psi = c*hpi*i_vals(2:end).^eps; % time integration
        
    case 4
        
        % cell death rate: time-dependent, input-dependent - d*k*hpi^(k-1)*i^eps,
        k = alpha(1);
        d = alpha(2);
        eps = alpha(3);
        int_psi = d*hpi^k*i_vals(2:end).^eps; % time integration
        
end

for n=1:length(moi)
    
    this_moi = moi(n);
    
    switch which_distribution
        
        case 1 % Poisson
            
            % probability cell not infected
            init_prob_moi0 = poisspdf(i_vals(1), this_moi);
            
            % probability of i virions entering - according to Poisson
            init_prob_moi1plus = poisspdf(i_vals(2:end), this_moi);
            
        case 2 % negative binomial
            
            matlab_r = k_dispersion;
            
            matlab_p = 1-this_moi/(matlab_r + this_moi); % formulation of negative binomial in matlab speak
            
            % probability cell not infected
            init_prob_moi0 = nbinpdf(i_vals(1),matlab_r,matlab_p);
            
            % probability of i virions entering - according to negative binomial
            init_prob_moi1plus = nbinpdf(i_vals(2:end),matlab_r,matlab_p);
            
    end
    
    init_prob_moi1plus(end) = init_prob_moi1plus(end) + (1 - (sum(init_prob_moi1plus)+init_prob_moi0));
    
    % int_psi: from integrating cell death rate model (above)
    prob_infectedcells_alive = sum(init_prob_moi1plus.*exp(-int_psi));
    
    % total prob alive
    prob_alive = init_prob_moi0 + prob_infectedcells_alive;
    
    % assume that cells only die if they are infected
    percent_alive(n) = 100*prob_alive;
    
end
