function u_best = pdMean(X,initial,M,N)
%% Compute the Frechet mean of all baseline matrices

%% algorithm parameters
alpha = 0.01; % updating step
IterNum = 100;
annealing_factor = 0.5; % shrink updating step

%% Computing Frechet mean
iter = 1;
u = initial;
u_best = u;
energy_best = 10^20;

while iter < IterNum
    ul = u;
    
    d = zeros(N,N);

    energy = 0;
    for i = 1:M
        dx = logMap(u,X{i});
    	d  = d + dx;
        % energy = energy + norm(dx)^2; % use an approximated Euclidean
        % metric
        energy = energy + pdDist(u,expMap(u,dx)); % use Riemannian metric
    end

    if (energy < energy_best)
        energy_best = energy;
        u_best = u;
    else
        alpha = alpha * annealing_factor;
        u = u_best;
        iter = iter + 1;
        %fprintf('%d\n',iter);
        continue; 
    end
    
    u = expMap(u,d * alpha);
    
    update = expMap(ul,u);
    
    if (mod(iter,20)==0)
        fprintf('Iteration %d update: %f, energy %f\n',iter, norm(update(:)), energy);
    end
    iter = iter + 1;
    
    %% optimal shrinking
% 	if (mod(iter,10) == 0)
%         alpha = alpha * 0.1;
% 	end
end

u_best = real(u_best);
end