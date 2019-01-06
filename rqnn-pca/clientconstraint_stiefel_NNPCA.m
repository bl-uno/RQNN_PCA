function [data,W] = clientconstraint_stiefel_NNPCA(D, k, methodoptions)
% Author=Bruno Losseau              Date=6/1/2019
% Writing it under the framework of:
% Optimization-on-manifolds-with-extra-constraints-master by Liu and Boumal
% available at https://github.com/losangle.
% INPUT. 
%    - D has to be non-negative, input data a matrix of size m x n
%    - k is a positive integer, number of clusters
%    - methodoptions
% OUTPUT.
%    - data, statistics containing
%               Row1: Maxviolation of orthogonality at solution
%               Row2: Cost, final solution 
%               Row3: Time to find the final solution
%   - W, output data, a matrix of size m x k


    data = NaN(3, 1); 
    [m,~] = size(D);
    M = stiefelfactory(m,k);

    problem.M = M;
    problem.cost = @(Y) costFun(Y,D);
    problem.egrad = @(Y) gradFun(Y,D);
    problem.D = D;
    [x0,~,~] = svds(D,k); % we must start with the svd solution, and do the
                          % flip to push the solution to be as non-neg as possible
for j=1:k
    negidx = x0(:,j)<0;
    isNegNormGreater = norm(x0(negidx,j),'fro') > norm(x0(~negidx,j),'fro');
    if isNegNormGreater
        x0(:,j) = -x0(:,j);
    end  
end

%   DEBUG only

%   checkgradient(problem);

%   -------------------------Set-up Constraints-----------------------
%   Nonnegativity of all entries

    ineq_constraints_cost = cell(1,1);
    ineq_constraints_cost{1} = @(U) 1/2*sum(sum(U(U<0).^2)); %Frobenius norm of the negative elements

    ineq_constraints_grad = cell(1,1);
    ineq_constraints_grad{1} = @(U)  min(U,0) ; 

    problem.ineq_constraint_cost = ineq_constraints_cost;
    problem.ineq_constraint_grad = ineq_constraints_grad;

    condet = constraintsdetail(problem);
    
%     Debug Only
%     checkconstraints(problem)

%     ------------------------- Solving ---------------------------
    options = methodoptions;
    
        %Quadratic
        fprintf('Starting Quadratic \n');
        timetic = tic();
        
        [W,~] = exactpenaltyViaQuadratic(problem, x0, D, options);
        time = toc(timetic);
        [maxviolation, meanviolation, cost] = evaluation(problem, W, condet);
        data(1) = max(maxviolation, manifoldViolation(W));
        data(2) = cost;
        data(3) = time;
    
        
        
        
     %------------------------sub functions-----------     
    function f = costFun(Y,D)
        f = -0.5*norm(Y'*D,'fro')^2; 
       %f = 0.5*norm(D-Y*(Y'*D),'fro')^2; %the same problem, under min and max forms
                                        %the second is slower per inner
                                        %iteration as we need to make more matrix manipulations
                                        % first is like PCA, second is like PNMF
   
    end

    function val = gradFun(Y,D)
        
        val = -D*(D'*Y);
    end

    function manvio = manifoldViolation(x)
        size(x)
        manvio = max(max(abs(x.'*x-eye(k))));
    end

        

end