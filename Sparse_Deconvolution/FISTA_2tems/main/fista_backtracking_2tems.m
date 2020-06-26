function [X1,X2,cost_iter] = fista_backtracking_2tems(calc_f, grad, Xinit1,Xinit2, opts, calc_F)   
% function [X, iter, min_cost] = fista_backtracking(calc_f, grad, Xinit, opts, calc_F)   
% * A Fast Iterative Shrinkage-Thresholding Algorithm for 
% Linear Inverse Problems: FISTA (backtracking version)
% * Solve the problem: X = arg min_X F(X) = f(X) + lambda*g(X) where:
%   - X: variable, can be a matrix.
%   - f(X): a smooth convex function with continuously differentiable 
%       with Lipschitz continuous gradient `L(f)` (Lipschitz constant of 
%       the gradient of `f`).
%  INPUT:
%       calc_f  : a function calculating f(x) in F(x) = f(x) + g(x) 
%       grad   : a function calculating gradient of f(X) given X.
%       Xinit  : a matrix -- initial guess.
%       opts   : a struct
%           opts.lambda  : a regularization parameter, can be either a scalar or
%                           a weighted matrix.
%           opts.max_iter: maximum iterations of the algorithm. 
%                           Default 300.
%           opts.tol     : a tolerance, the algorithm will stop if difference 
%                           between two successive X is smaller than this value. 
%                           Default 1e-8.
%           opts.verbose : showing F(X) after each iteration or not. 
%                           Default false. 
%           opts.L0 : a positive scalar. 
%           opts.eta: (must be > 1). eta in the algorithm (page 194)

%       calc_F: optional, a function calculating value of F at X 
%               via feval(calc_F, X). 
%  OUTPUT:
%      X        : solution
% ************************************************************************
% * Date created    : 11/06/17
% * Author          : Tiep Vu 
% * Date modified   : 
% ************************************************************************

    if ~isfield(opts, 'max_iter')
        opts.max_iter = 500;
    end
    if ~isfield(opts, 'regul')
        opts.regul = 'l1';
    end     
    if ~isfield(opts, 'pos')
        opts.pos = false;
    end
    
    if ~isfield(opts, 'tol')
        opts.tol = 1e-8;
    end
    
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end

    if ~isfield(opts, 'L0')
        opts.L0 = 1;
    end 
    if ~isfield(opts, 'eta')
        opts.eta = 2;
    end 
%     if ~isfield(opts, 'lambda2_scaler')
%         opts.lambda2_scaler = 4;
%     end 

    %% projection 
    function res = projection(U, opts) 
        switch opts.regul
            case 'l1'
                res = proj_l1(U, opts);
            case 'l2' 
                res = proj_l2(U, opts);
            case 'l12'
                res = proj_l2(U', opts)';
        end
    end 
    %% computer g 
    function res = g(x1,x2) 
        switch opts.regul
            case 'l1'
                res = sum(sum(abs(opts.lambda1.*x1)))+sum(sum(abs(opts.lambda2.*x2)));
                %res = opts.lambda1.*norm1(x1)+opts.lambda2.*opts.lambda2_scaler.*norm1(x2)-opts.lambda2.*(opts.lambda2_scaler-1).*norm1(x2.*(x1~=0));
%             case 'l2' 
%                 res = opts.lambda*norm2_cols(x);
%             case 'l12'
%                 res = opts.lambda*norm12(x);
        end
    end 

    %% computer Q 
    function res = calc_Q(x1, x2, y1, y2, L,g1,g2) 
        % based on equation 2.5, page 189
        res = feval(calc_f, y1, y2) + (x1 - y1)'*g1 +(x2 - y2)'*g2 ...
                    + L/2*normF2(x1 - y1) + L/2*normF2(x2 - y2) + g(x1,x2);
    end 

    x1_old = Xinit1;
    x2_old = Xinit2;
    y1_old = Xinit1;
    y2_old = Xinit2;
    t_old = 1;
    iter = 0;
    cost_old = 1e10;
    %% MAIN LOOP
    opts_proj = opts;
%     opts_proj.lambda = lambdaLiv;
    L = opts.L0;
    opts0 = opts;
    cost_iter=zeros(opts.max_iter,1);
    if parallel.gpu.GPUDevice.isAvailable
        cost_iter=gpuArray(cost_iter);
    end
    while  iter < opts.max_iter
        if opts.verbose
            tic
        end
        iter = iter + 1;
        % find i_k 
        Lbar = L; 
        [g1,g2]=feval(grad, y1_old, y2_old);
        while true 
            opts0.lambda = opts.lambda1/Lbar; 
            zk1 = projection(y1_old - 1/Lbar*g1, opts0);
            opts0.lambda = opts.lambda2/Lbar;
            zk2 = projection(y2_old - 1/Lbar*g2, opts0);
            F = feval(calc_F, zk1, zk2);
            Q = calc_Q(zk1, zk2, y1_old, y2_old, Lbar,g1,g2);
            if F <= Q 
                break;
            end
            Lbar = Lbar*opts.eta; 
            L = Lbar; 
        end 
%         L = L/opts.eta;
        opts_proj.lambda = opts.lambda1/L;
        x1_new = projection(y1_old - 1/L*g1, opts_proj);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y1_new = x1_new + (t_old - 1)/t_new * (x1_new - x1_old);
    
        opts_proj.lambda = opts.lambda2/L;
        x2_new = projection(y2_old - 1/L*g2, opts_proj);
        y2_new = x2_new + (t_old - 1)/t_new * (x2_new - x2_old);
        %% check stop criteria
        e = (norm1(x1_new - x1_old)+norm1(x2_new - x2_old))/(numel(x1_new)+numel(x2_new));
        if e < opts.tol
            break;
        end
        %% update
        t_old = t_new;
        x1_old = x1_new;
        y1_old = y1_new;
        x2_old = x2_new;
        y2_old = y2_new;
        cost_new = feval(calc_F, x1_new, x2_new);
        cost_iter(iter)=cost_new;
        %% show progress
        if opts.verbose
            if nargin ~= 0
                if cost_new <= cost_old 
                    stt = 'YES.';
                else 
                    stt = 'NO, check your code.';
                end
                fprintf('iter = %3d, cost = %f, cost decreases? %s\n', ...
                    iter, cost_new, stt);
                cost_old = cost_new;
                toc
            else 
                if mod(iter, 5) == 0
                    fprintf('.');
                end
                if mod(iter, 10) == 0 
                   fprintf('%d', iter);
                end     
            end        
        end 
    end
    X1 = x1_new;
    X2 = x2_new;
end 