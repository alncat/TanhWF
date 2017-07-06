%%  Compute the truncated gradient based on the Poisson log-likelihood function

function [grad, corr, var] = compute_grad(z, y, Params, A, At, normest, t)
    m = length(y);
    corr = 0;
    var = 0;
    yz = A(z);
    Kt = 1/m* norm(abs(yz(:)).^2 - y(:), 1); 
    
    if strcmp(Params.grad_type,'TWF_Poiss') == 1   % truncated gradient / Wirtinger flow
        % truncation rules
        Eub =  abs(yz) / norm(z)   <= Params.alpha_ub;
        Elb =  abs(yz) / norm(z)   >= Params.alpha_lb;
        Eh  =  abs(y - abs(yz).^2) <= Params.alpha_h * Kt / norm(z) * abs(yz);
        
        grad  = 1/m* At( 2* ( abs(yz).^2-y ) ./ (abs(yz).^2) .*yz ...
                          .* Eub .* Elb .* Eh );    % truncated Poisson gradient
                      
    elseif strcmp(Params.grad_type,'WF_Poiss') == 1    % untruncated gradient / Wirtinger flow
        grad  = 1/m* At( 2* ( abs(yz).^2-y ) ./ (abs(yz).^2) .*yz ); % Poisson gradient
    elseif strcmp(Params.grad_type, 'LIN') == 1
        %norm_y = sqrt(y) / normest;
        %norm_z = yz / norm(z);
        var = (yz.^2 + y - 2*abs(yz .* sqrt(y)));
        %w = max([tanh(sqrt(y) ./ abs(yz))';tanh(abs(yz) ./ sqrt(y))']);
        %w = w' .* sign(yz);
        w = tanh((yz .* sqrt(y))./abs(var));
        %w = sign(yz);
        grad = 2/m* At( yz - w .* sqrt(y));
    elseif strcmp(Params.grad_type,'MX_LLK') == 1 % cosh type likelihood
        %var = sum((y + (yz).^2 - 2*tanh((yz .* sqrt(y)) ./ abs(var)).*(yz .* sqrt(y))))/length(y);
        if t < 25
            %var = norm(y - abs(yz).^2, 1)/length(y);
            var_d = [sqrt(y) - yz; sqrt(y) + yz];
            var_d = var_d.^2;
            var = sum(var_d/2)/length(y);
            %disp([var, t])
            % truncation rules
%             Eub =  abs(yz) / norm(z)   <= Params.alpha_ub;
%             Elb =  abs(yz) / norm(z)   >= Params.alpha_lb;
%             Eh  =  abs(wy - abs(yz).^2) <= Params.alpha_h * Kt / norm(z) * abs(yz);
%             tmp  = 1/m* At( 2* ( abs(yz).^2-y ) ./ (abs(yz).^2) .*yz ...
%                           .* Eub .* Elb .* Eh );    % truncated Poisson gradient
%             var = max(0.01*sum(tmp.^2), 30);
        else
            var_d = [(sqrt(y) - yz)'; (sqrt(y) + yz)'];
            var_d = var_d.^2;
            var = min(var_d);
            %if t < 50
            var = sum(var)/length(y)*8;
            %else
            %var = norm(z)^2 + normest^2 - 2*sum(abs(yz).*sqrt(y))/length(y);
            %var = var*4;
            %end
        end
        corr = sum(yz .* sqrt(y))/length(y);
        w = tanh((yz .* sqrt(y)) ./ abs(var));
        %sum(w)/length(y);
        grad = 2/m* At( yz - w .* sqrt(y));
        a = max(abs(grad));
    end