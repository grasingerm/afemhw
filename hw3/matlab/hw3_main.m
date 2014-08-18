# program constants

MAX_ITERATIONS = 1e5;
DELTA_LAMBDA = 0.05;
TOLERANCE = 1e-3;

#{
    Read initial stresses and strains
    Define initial guess u_o, lambda_o
#}

u_o = zeros(8,1);
lambda = 0;
load_cycles = 1/DELTA_LAMBDA;
total_iterations = 0;

disp("Maximum iterations: "), disp(MAX_ITERATIONS)
disp("Load cycles: "), disp(load_cycles)
disp("Error tolerance: "), disp(TOLERANCE)

for i=1:load_cycles
    lambda = lambda + DELTA_LAMBDA
    k = 0;
    
    while k < MAX_ITERATIONS && error > TOLERANCE
        k++;
        total_iterations++;
        
        disp("total iterations "), disp(total_iterations)
        disp("local iteration "), disp(k)
        
        # compute residual and jacobian matrix
        
        
        # enforce homogenous boundary conditions
        
        # compute correction
        
        # update solution
        
        # compute error
        
    end
    
    # make updates
    
end
