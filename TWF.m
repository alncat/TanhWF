%% Implementation of the truncated Wirtinger Flow (TWF) algorithm proposed in the paper
%  ``Solving Random Quadratic Systems of Equations Is Nearly as Easy as Solving Linear Systems'' by Y. Chen and E. J. Candès.
%  The code below is adapted from implementation of the Wirtinger Flow algorithm designed and implemented by E. Candes, X. Li, and M. Soltanolkotabi

function [Relerrs, corrs, vars] = TWF(y, x, Params, A, At)    
%% Initialization
    npower_iter = Params.npower_iter;           % Number of power iterations 
    z0 = randn(Params.n1,Params.n2); z0 = z0/norm(z0,'fro');    % Initial guess 
    normest = sqrt(sum(y(:))/numel(y(:)));    % Estimate norm to scale eigenvector  
    
    for tt = 1:npower_iter,                     % Truncated power iterations
        yz = A(z0);
        %varest = sum(y - abs(yz).^2)/length(y);
        if strcmp(Params.init, 'Trunc') == 1
            ytr = y.* (abs(y) <= Params.alpha_y^2 * normest^2 );
        else
            ytr = tanh(y /(4 * normest^2));
            %ytr = 1 - exp(-y/(2 * normest^2));
        end
        z0 = At( ytr .* yz ); z0 = z0/norm(z0,'fro');
    end
    %z = x + 1.5*randn(Params.n1, Params.n2);
    %z = normest*z/norm(z,'fro');
    z = normest * z0;                   % Apply scaling 
    Relerrs = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro'); % Initial rel. error

    %% Loop
    grad_type = Params.grad_type;
    if strcmp(grad_type, 'TWF_Poiss') == 1
        mu = @(t) Params.mu; % Schedule for step size 
    elseif strcmp(grad_type, 'WF_Poiss') == 1
        tau0 = 330;                         % Time constant for step size
        mu = @(t) min(1-exp(-t/tau0), 0.2); % Schedule for step size  
    elseif strcmp(grad_type, 'MX_LLK') == 1
        tau0 = 330;
        mu = @(t) min(1-exp(-t/tau0), 0.9);
        %mu = @(t) Params.mu;
    end
    
    var = Params.alpha_y*normest^2;
    v = zeros(size(z));
    momentum = 0.7;
    l_r = 1.5e-1;
    vars = [];
    corrs = [];
    diffs = [];
    for t = 1: Params.T,
        %var = Relerrs(t)*norm(x, 'fro');
        v_pre = v;
        %var_pre = var;
        [grad, corr, var] = compute_grad(z, y, Params, A, At, normest, t);
        corrs = [corrs, abs(z'*x)/(norm(z)*norm(x))];
        vars = [vars, var];
        v = momentum*v - l_r * grad;
        z_pre = z;
        z = z - momentum*v_pre + (1+momentum)*v;
        diffs = [diffs, norm(z - z_pre,'fro')];
        %var = sum(grad.^2);
        %z = z - mu(t) * grad;             % Gradient update 
        Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')];
        %Relerrs = [Relerrs, min(norm(x - z)/norm(x, 'fro'), norm(x+z)/norm(x, 'fro'))];
    end
