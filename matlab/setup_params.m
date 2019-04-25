function params = setup_params(theta,params)
%SETUP_PARAMS Assign parameter vector to PARAMS struct or get parameter bounds.

% Also return bounds for fitting
if isempty(theta)
    bounds.LB = [];
    bounds.UB = []; 
    bounds.PLB = [];
    bounds.PUB = []; 
    bounds.x0 = [];     
end

iParam = 1;
while iParam <= numel(params.names)
    pname = params.names{iParam};
    if ~isfield(params,pname) || isempty(params.(pname))
        error(['Parameter ''' pname ''' not in PARAMS struct.']);
    end
    pnum = numel(params.(pname));
    % Check that vector parameters are passed correctly
    if pnum > 1
        for jParam = iParam+1:iParam+pnum-1
            if ~strcmp(params.names{jParam},pname)
                error(['Parameter ''' pname ''' should span multiple elements.']);
            end
        end
    end
    
    if ~isempty(theta)
        params.(pname)(1:pnum) = theta(iParam:iParam+pnum-1);

        % Extra parameter processing
        switch pname
            case 'sigma'
                Nsig = numel(params.sigma);
                params.sigma = exp(params.sigma);
                if Nsig == numel(params.sigma_contrasts)
                    params.sigma_poly = polyfit(log(params.sigma_contrasts),params.sigma,numel(params.sigma_contrasts)-1);
                elseif Nsig == 2*numel(params.sigma_contrasts)
                    params.sigma_poly_left = polyfit(log(params.sigma_contrasts),params.sigma(1:Nsig/2),numel(params.sigma_contrasts)-1);                
                    params.sigma_poly_right = polyfit(log(params.sigma_contrasts),params.sigma(Nsig/2+1:end),numel(params.sigma_contrasts)-1);                                    
                end
            case 'precision'
                params.precision = exp(params.precision);                
            case 'attention_factor'
                params.attention_factor = exp(params.attention_factor);
            case 'softmax_eta'
                params.softmax_eta = exp(params.softmax_eta);
            case 'runlength_tau'
                params.runlength_tau = exp(params.runlength_tau);
                params.runlength_prior = ['@(t) exp(-t/' num2str(params.runlength_tau,'%.8f') ')'];
            case 'runlength_min'
                params.runlength_min = exp(params.runlength_min);
            case 'psycho_sigma'
                params.psycho_sigma = exp(params.psycho_sigma);
            case 'prob_low'
                % Assign low probability to probability vector
                [~,idx_low] = min(params.p_true_vec);
                params.p_vec(idx_low) = params.prob_low;
                % Change high probability symmetrically 
                % (if PROB_HIGH is not set independently)
                params.prob_high = 1 - params.prob_low;
                [~,idx_max] = max(params.p_true_vec);
                params.p_vec(idx_max) = 1 - params.prob_low;
            case 'prob_high'
                [~,idx_max] = max(params.p_true_vec);
                params.p_vec(idx_max) = params.prob_high;
        end
    end
    
    % Return parameter bounds
    if isempty(theta)
        % BVEC is LB,UB,PLB,PUB,X0 for that parameter
        switch pname
            case 'sigma'
                bvec = [log(0.1),log(180),log(2),log(60),log(10)];
            case 'precision'
                bvec = [log(1/180^2),log(1),log(1/40^2),log(1/5^2),log(1/15^2)];
            case 'precision_power'
                bvec = [0,10,0.1,0.9,0.5];                
            case 'lapse_rate'
                bvec = [0,1,0.01,0.2,0.05];
            case 'lapse_bias'
                bvec = [0.01,0.99,0.1,0.9,0.5];
            case 'attention_factor'
                bvec = [-2,2,-1,1,0];                
            case 'softmax_eta'
                bvec = [-20,20,-10,10,0];
            case 'softmax_bias'
                bvec = [-50,50,-10,10,0];
            case 'runlength_tau'
                bvec = [0,log(200),log(10),log(100),log(60)];
            case 'runlength_min'
                bvec = [0,log(100),log(5),log(40),log(20)];                
            % Psychometric function parameters
            case 'psycho_mu'
                bvec = [-1,1,-0.2,0.2,0];
            case 'psycho_sigma'
                bvec = [log(0.001),log(10),log(0.01),log(1),log(0.1)];
            case {'psycho_gammalo','psycho_gammahi'}
                bvec = [0,0.5,0.01,0.1,0.05];
            case {'prob_low'}
                bvec = [1e-6,0.5,0.05,0.45,0.2];
            case {'prob_high'}
                bvec = [0.5,1-1-e6,0.55,0.95,0.8];
        end
        
        bounds.LB = [bounds.LB,bvec(1)*ones(1,pnum)];
        bounds.UB = [bounds.UB,bvec(2)*ones(1,pnum)]; 
        bounds.PLB = [bounds.PLB,bvec(3)*ones(1,pnum)];
        bounds.PUB = [bounds.PUB,bvec(4)*ones(1,pnum)]; 
        bounds.x0 = [bounds.x0,bvec(5)*ones(1,pnum)]; 
    end
    
    iParam = iParam + pnum;
end

if isempty(theta)
    params = bounds;
else
    params.theta = theta;
end

end