function [xfinal,info] = exactpenaltyViaQuadratic(problem0, x0, D, options)
    % Writing it under the framework of
    % Optimization-on-manifolds-with-extra-constraints-master by Liu and Boumal
    % available at https://github.com/losangle.

    condet = constraintsdetail(problem0);
    info = 0;
   %Standard Setting, you should only check these parameters!
    localdefaults.rho = 1e2; 
    localdefaults.C =1.1; %between 1.005 and 1.5
    localdefaults.beta = 1.00; %between 1 and 1.2, usually no more than 1.01
    localdefaults.tau = sqrt(sum(sum(D.*D)));
    localdefaults.relGradConvTreshold = 1e2; % relative to D
    localdefaults.relNegativityTreshold = 1e-4; 
    localdefaults.maxOuterIter = 2e3;
    localdefaults.maxInnerIter = 2e3; %there shouldn t be more than 2 iters on average, exceptionaly 20 
    localdefaults.verbosity = 0;
    localdefaults.showErr= 0;
    
    
    %Inner Loop Setting
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    
    rho = options.rho;
    C = options.C;
    beta = options.beta;
    tau = options.tau;
    relNegativityTreshold = options.relNegativityTreshold;
    relGradConvTreshold = options.relGradConvTreshold;
    M = problem0.M;
    xCur = x0;
    showErr = options.showErr;
    
    errfit = zeros(options.maxOuterIter,1);
    errneg = zeros(options.maxOuterIter,1);
    totaltime = tic();
    
    % Outer loop
    for OuterIter = 1 : options.maxOuterIter
        costfun = @(X) cost_exactpenalty(X, problem0, rho);
        gradfun = @(X) grad_exactpenalty(X, problem0, rho);
        problem.cost = costfun;
        problem.grad = gradfun;
        problem.M = M;
        inneroptions.tolgradnorm = tau*relGradConvTreshold/beta;
        inneroptions.verbosity = options.verbosity;
        inneroptions.maxiter = options.maxInnerIter;
        inneroptions.minstepsize = options.minstepsize;
       % inneroptions.linesearch = @linesearch; 
   
    
        % Inner loop using steepest descent
        [xCur] = steepestdescent(problem, xCur, inneroptions);
        errfit(OuterIter) = computeFitError(D,xCur);
        relNorm= norm(xCur(xCur<0),'fro')/norm(xCur,'fro');
        errneg(OuterIter)=relNorm;        
        
        % Updating the penalty parameter
        rho = rho*C;
         
         
        if mod(OuterIter,10)==0
            fprintf('Iteration: %d    ', OuterIter); 
            fprintf('\n Neg: %.3e, misft %.3e, rho %.3e \n', relNorm, errfit(OuterIter),rho)
        end

        
        % Check if we have reached sufficient negativity
        if    relNorm<=relNegativityTreshold % Focus is only on non-negativity, we could add a parameter to check for 
                                             % convergence to a local infeasible solution
            break;
        end
        
         
        
        if toc(totaltime) > options.maxtime
            break;
        end
        
    end

    if showErr==1 %show evolution of the error on negativity and fit
        hf = figure(1); hold on;hp = semilogy(errneg); grid on; axis square; set(hf,'Color','white'); set(hp,'LineWidth',2.0)
        hf = figure(2); hold on; hp = semilogy(errfit); grid on; axis square; set(hf,'Color','white'); set(hp,'LineWidth',2.0)
    end
        
    xfinal = xCur;

    function val = cost_exactpenalty(x, problem, rho)
        val = getCost(problem, x);
        % Adding ineq constraint cost
        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(x);
                val = val + rho * cost_at_x;
            end
        end
    end

    function val = grad_exactpenalty(x, problem, rho)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                if (cost_numineq > 0)
                    gradhandle = problem.ineq_constraint_grad{numineq};
                    constraint_grad = gradhandle(x);
                    constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                    val = problem.M.lincomb(x, 1, val, rho, constraint_grad);
                end
            end
        end
     
    end
end

function fitErr = computeFitError(D,U)
    UtD = U'*D;
    fitErr =abs(sum(sum(D.*D))-sum(sum(UtD.*UtD)))/sum(sum(D.*D));
end

