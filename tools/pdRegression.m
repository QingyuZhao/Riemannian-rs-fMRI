function [p_best,v_best] = pdRegression(X)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% implementation of geodesic regression on positive define cone based on
% approximated gradient descent described in the following
%  
% Hyunwoo J Kim, et al. Multivariate general linear models (MGLM) on
% Riemannian manifolds with applications to statistical analysis of
% diffusion weighted images, CVPR 2014

    %% regression algorithm parameters
    IterNum = 101;
    alpha = 0.1; % update step for p (starting position)
    beta = 0.1; % update step for v (shooting vector)
    annealing_alpha = 0.5; % shrink update steps (alpha, beta)
    annealing_beta= 0.25; % shrink update steps (alpha, beta)
    
    convergence_thre_p = 0.00001;
    convergence_thre_v = 0.001;
    
    %% initialization
    M = length(X);
    N = size(X{1},1);
    
    iter = 1;
    p = X{1}; % initialize baseline point
    v = zeros(N,N);
    for i = 2:M
        v = v + logMap(X{1},X{i})/(i-1);
    end
    v = v / M; % initialize shooting direction
    
    p_best = p;
    v_best = v;
    energy_best = 10^20;
    
    while iter < IterNum
        p_last = p;
        v_last = v;
        
        %% compute regression energy
        energy = 0;
        update_p = zeros(N,N);
        for i = 1:M
            regressed_p{i} = expMap(p,v*(i-1));
            Jacobi{i} = logMap(regressed_p{i},X{i});
            Jacobi_base{i} = pdParaTrans(Jacobi{i},regressed_p{i},p);
            
            update_p = update_p + Jacobi_base{i};
            energy = energy + pdDist(regressed_p{i},X{i});
        end
        
        %% update best parameters
        if (mod(iter,50) == 0)
            fprintf('Iter %d: Energy: %f \n',iter,energy);
        end
        if (energy < energy_best)
            energy_best = energy;
            p_best = p;
            v_best = v;
        else
            p = p_best;
            v = v_best;
            alpha = alpha * annealing_alpha;
            beta = beta * annealing_beta;
            iter = iter + 1;
            continue; 
        end
        
        %% update p
        p = expMap(p,update_p * alpha);

        %% update v
        update_v = zeros(N,N);
        for i = 2:M
            regressed_p{i} = expMap(p,v*(i-1));
            Jacobi{i} = logMap(regressed_p{i},X{i});
            Jacobi_base{i} = pdParaTrans(Jacobi{i},regressed_p{i},p);
            
            update_v = update_v + Jacobi_base{i}*(i-1);
            energy = energy + pdDist(regressed_p{i},X{i});
        end

        v = v + update_v * beta;
        
        %% stop criteria when updates are small enough
        error1 = logMap(p_last,p);
        error2 = v_last - v;
        
        error1 = norm(error1(:));
        error2 = norm(error2(:));
        
        if (mod(iter,10) == 0)
            fprintf('       iter %d error1 %f    error2 %f\n',iter, norm(error1(:)),norm(error2(:)));
        end
        
        if (error1 < convergence_thre_p) && (error2 < convergence_thre_v)
            break;
        end
        
        %% optional annealing to help more stable (but slower convergence) on specific dataset
%         if (mod(iter,20) == 0) && (iter <= 60)
%             alpha = alpha * 0.5;
%             beta = beta * 0.25;
%         end
        
        iter = iter + 1;
    end
end