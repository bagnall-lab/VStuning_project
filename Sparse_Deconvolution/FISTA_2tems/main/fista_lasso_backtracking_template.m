function X = fista_lasso_backtracking_template(Y, template, Xinit, opts)
%% Y is n*1 matrix,template is m*1 matrix, therefore X is (n-m+1)*1
    if ~isfield(opts, 'backtracking')
        opts.backtracking = false;
    end 
    opts.regul = 'l1';
    opts = initOpts(opts);
    lambda = opts.lambda;

%     if numel(lambda) > 1 && size(lambda, 2)  == 1
%         lambda = repmat(opts.lambda, 1, size(Y, 2));
%     end
    if numel(Xinit) == 0
        Xinit = zeros(size(Y,1)-size(template,1)+1, size(Y,2));
    end
    %% cost f
    function cost = calc_f(X)
         cost = 1/2 *normF2(Y - conv(X,template));
    end 
    %% cost function 
    function cost = calc_F(X)
        if numel(lambda) == 1 % scalar 
            cost = calc_f(X) + lambda*norm1(X);
         
        elseif numel(lambda) == numel(X)
            cost = calc_f(X) + norm1(lambda.*X);
        elseif numel(lambda) == size(X, 1) 
            lambda1 = repmat(lambda, 1, size(size(X, 2)));
            cost = calc_f(X) + norm1(lambda1.*X);
        end
    end 
    %% gradient
    %DtD = D'*D;
    %DtY = D'*Y;
    function res = grad(X) 
        %res = DtD*X - DtY(:, i);
        E=conv(X,template)-Y;
        TemE_corr=xcorr(template,E);
        res=TemE_corr(length(E):-1:length(template));
%         res=zeros(size(X,1),size(X,2));
%         for k=1:size(res,2)
%             for j=1:size(res,1)
%                 res(j,k)=template'*E(j:j+size(template,1)-1,k);
%             end
%         end
    end 
    % Checking gradient 
    if opts.check_grad
        check_grad(@calc_f, @grad, Xinit);
    end 

    opts.max_iter = 500;
    % for backtracking, we need to optimize one by one 
    X = zeros(size(Xinit));
    for i = 1:size(X, 2) 
       X(:, i) = fista_backtracking_template(@calc_f, @grad, Xinit(:, i), opts, ...
                                        @calc_F);
    end 
end 