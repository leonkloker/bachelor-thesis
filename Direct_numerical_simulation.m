% Script that runs a direct numerical simulation of Darcy's equation (10), the
% continuity equation (9) and the convection-diffusion equation of salt (11) with
% the corresponding boundary conditions (13). 
% 
% The dimensionless height alpha, wavenumber a, Rayleigh number Ra,
% amount of vertical and lateral cells N and M, respectively, end time, 
% timestep and perturbation seed can be adjusted.
%
% The calculated state vector at time t = (i-1)*dt is saved in evolution(i,:).
% The state can be visualized either by using the function plot_dns or by
% using the transform function to get a struct containing the 2d arrays
% 'c', 'p', 'vx' and 'vz'. These arrays can then be plotted using e.g. surf or
% pcolor.
%
% The development of the perturbation amplitude \sigma_\chi is saved in instability and
% ban be plotted using semilogy(instability).

% Dimensionless height
alpha = 5;

% Wavenumber of the perturbation mode to be investigated
a = 0.26;

% Rayleigh number
Ra = 14;

% Cells in vertical direction (the amount of vertical cells should be increased
% for large alpha; reference: N = 360 for alpha = 5 is needed to accurately resolve 
% the perturbations)
N = 100;

% Cells in lateral direction (10 are already sufficient as only one
% wavelength of the perturbation has to be resolved)
M = 8;

% End time of the simulation
T = 5;
% Timestep
dt = 0.01;

% Vertical size of finite volume cells (equisized by default, but can be 
% changed in order to resolve the top of the porous medium more finely)
hz = ones(N+2,1)*alpha/N;

% Corresponding wavelength / domain width
lambda = 2*pi/a;
% Lateral size of the finite volume cells (only equisized cells possible)
hx = lambda/M;

% Initialize time to zero
t = 0

% Initialize the perturbation seed of the salt concentration
u = initialize_state_cosine(M, N, alpha, a, hz, 10^-10);

% Alternative: u = initialize_state_custom(M, N, alpha, Ra, a, tc, hz, amplitude)
% This initialization only works if tc is the time of onset for the given 
% alpha, Ra and a, otherwise the Chebyshev-Galerkin method does not yield 
% a non-zero perturbation

% Array saving the standard deviation \sigma_\chi at each timestep
instability = zeros(cast(T/dt,'int32')+1,1);
u1 = transform(M,N,u);
instability(1) = std(u1.c(2,2:end-1));

% Array saving the entire state vector at each timestep
evolution = zeros(cast(T/dt,'int32')+1,max(size(u)));
evolution(1,:) = u;

% Integrate the equations
while t + dt < T + 10^-9

    % Get the coefficient matrix of the implicit Euler scheme
    A = equation_system(M, N, Ra, dt, hx, hz, u);
    % Get the load vector
    b = load(M, N, Ra, u);

    % Save the old state
    u_old = u;
    % Solve the equation system
    u = A\b;

    % Iteration counter for the nonlinear flux term
    i = 1;
    
    % Iterate until the state difference is smaller than the threshold
    while norm(u_old-u)/max(size(u)) > 10^-9 && i < 100
        i = i+1;
        A = update_equation_system(A, M, N, Ra, dt, hx, hz, u);
        u_old = u;
        u = A\b;
    end

    if i == 100
        error(sprintf('Iteration of advective flux term did not converge in %d steps!',i));
        return;
    end

    % Update time
    t = t + dt
    
    % Save \sigma_\chi and the current state 
    u1 = transform(M,N,u);
    perturbation_amplitude = std(u1.c(2,2:end-1))
    instability(cast(t/dt,'int32')+1) = perturbation_amplitude;
    evolution(cast(t/dt,'int32')+1,:) = u;

    % Add a new perturbation seed when the first one has become too small
    if perturbation_amplitude < 10^-14
        u = add_cosine(u, M, N, alpha, a, hz, 10^-8);
    end
end

function A = update_equation_system(A, M, N, Ra, dt, hx, hz, u_old)
% Updates the equation system during iteration of the advective flux term. 
% Only the salt concentration appearing in the nonlinear advection term is
% changed
%
% - A: matrix of the linear equation system to be solved at the previous iteration
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - Ra: Rayleigh number
% - dt: timestep
% - hx: lateral cell size
% - hz: vector containing vertical cell sizes
% - u_old: state vector after the previous iteration

    omega = 1;
    offset = 2*N*M+(N-1)*M;
    for i=1:N
        for j=1:M
            % Donor Cell scheme for the advection term
            A(offset+(i-1)*M+j,index('vx', i+1, j+1, M, N)) = -Ra*dt/hx * (u_old(index('c', i+1, j, M, N)) + u_old(index('c', i+1, j+1, M, N)) + omega * sign(u_old(index('vx', i+1, j+1, M, N))) * (u_old(index('c', i+1, j, M, N)) - u_old(index('c', i+1, j+1, M, N))))/2; 
            A(offset+(i-1)*M+j,index('vx', i+1, j+2, M, N)) = Ra*dt/hx * (u_old(index('c', i+1, j+1, M, N)) + u_old(index('c', i+1, j+2, M, N)) + omega * sign(u_old(index('vx', i+1, j+2, M, N))) * (u_old(index('c', i+1, j+1, M, N)) - u_old(index('c', i+1, j+2, M, N))))/2;
            A(offset+(i-1)*M+j,index('vz', i+1, j+1, M, N)) = Ra*dt/hz(i+1) * (u_old(index('c', i+1, j+1, M, N)) + u_old(index('c', i, j+1, M, N)) + omega * sign(u_old(index('vz', i+1, j+1, M, N))) * (u_old(index('c', i+1, j+1, M, N)) - u_old(index('c', i, j+1, M, N))))/2;
            A(offset+(i-1)*M+j,index('vz', i+2, j+1, M, N)) = -Ra*dt/hz(i+1) * (u_old(index('c', i+2, j+1, M, N)) + u_old(index('c', i+1, j+2, M, N)) + omega * sign(u_old(index('vz', i+2, j+1, M, N))) * (u_old(index('c', i+2, j+1, M, N)) - u_old(index('c', i+1, j+2, M, N))))/2;
        end
    end   
end

function b = load(M, N, Ra, u_old)
% Get the right-hand side of the linear equation system
%
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - Ra: Rayleigh number
% - u_old: state vector at previous timestep

    b = zeros(max(size(u_old)),1);
    offset = 2*N*M;
    offset = offset + M*(N-1);
    for i=1:N
        for j=1:M
            b(offset+(i-1)*M+j) = u_old(index('c', i+1, j+1, M, N));
        end
    end
    offset = offset + N*M;
    b(offset+1:offset+M) = 0;
    b(offset+M+1:offset+2*M) = 1/Ra;
    b(offset+2*M+1:offset+3*M) = -1;
    b(offset+3*M+1:offset+4*M) = 1/Ra;
end

function A = equation_system(M, N, Ra, dt, hx, hz, u_old)
% Get the matrix describing the linear equation system, that has to be
% solved at each timestep
%
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - Ra: Rayleigh number
% - dt: timestep
% - hx: lateral cell size
% - hz: vector containing vertical cell sizes
% - u_old: state vector at the previous timestep

    omega = 1;
    A = zeros(max(size(u_old)));

    for i=1:N
        for j=1:M
            if i == N && j == 1
                % Fix the pressure in one cell, such that equation system
                % is uniquely solvable
                A((i-1)*M+j,index('p', N+1, 1, M, N)) = 1;
            else
                % Continuity equation in all cells
                A((i-1)*M+j,index('vx', i+1, j+2, M, N)) = hz(i+1);
                A((i-1)*M+j,index('vx', i+1, j+1, M, N)) = -hz(i+1);
                A((i-1)*M+j,index('vz', i+2, j+1, M, N)) = -hx;
                A((i-1)*M+j,index('vz', i+1, j+1, M, N)) = hx;
            end
        end
    end

    offset = N*M;
    for i=1:N
        for j=1:M
            % Darcy's law in x-direction on interior edges
            A(offset+(i-1)*M+j, index('vx', i+1, j+2, M, N)) = 1;
            A(offset+(i-1)*M+j, index('p', i+1, j+2, M, N)) = 1/hx;
            A(offset+(i-1)*M+j, index('p', i+1, j+1, M, N)) = -1/hx;
        end
    end

    offset = offset + N*M;
    for i=1:N-1
        for j=1:M
            % Darcy's law in z-direction on interior edges
            A(offset+(i-1)*M+j, index('vz', i+2, j+1, M, N)) = 1;
            A(offset+(i-1)*M+j, index('c', i+1, j+1, M, N)) = hz(i+2)/(hz(i+1)+hz(i+2));
            A(offset+(i-1)*M+j, index('c', i+2, j+1, M, N)) = hz(i+1)/(hz(N+1)+hz(N+2));
            A(offset+(i-1)*M+j, index('p', i+1, j+1, M, N)) = 2/(hz(i+1)+hz(i+2));
            A(offset+(i-1)*M+j, index('p', i+2, j+1, M, N)) = -2/(hz(i+1)+hz(i+2));
        end
    end

    offset = offset + (N-1)*M;
    for i=1:N
        for j=1:M
            % Salt convection-diffusion equation in all cells
            A(offset+(i-1)*M+j,index('c', i+1, j+1, M, N)) = 1 + 2*dt*(1/(hz(i+1)*(hz(i)+hz(i+1))) + 1/(hz(i+1)*(hz(i+1)+hz(i+2)))) + 2*dt/(hx^2);
            A(offset+(i-1)*M+j,index('c', i, j+1, M, N)) = -2*dt/(hz(i+1)*(hz(i)+hz(i+1)));
            A(offset+(i-1)*M+j,index('c', i+2, j+1, M, N)) = -2*dt/(hz(i+1)*(hz(i+1)+hz(i+2)));
            A(offset+(i-1)*M+j,index('c', i+1, j, M, N)) = -dt/(hx^2);
            A(offset+(i-1)*M+j,index('c', i+1, j+2, M, N)) = -dt/(hx^2);

            % Donor Cell scheme for the advection term
            A(offset+(i-1)*M+j,index('vx', i+1, j+1, M, N)) = -Ra*dt/hx * (u_old(index('c', i+1, j, M, N)) + u_old(index('c', i+1, j+1, M, N)) + omega * sign(u_old(index('vx', i+1, j+1, M, N))) * (u_old(index('c', i+1, j, M, N)) - u_old(index('c', i+1, j+1, M, N))))/2; 
            A(offset+(i-1)*M+j,index('vx', i+1, j+2, M, N)) = Ra*dt/hx * (u_old(index('c', i+1, j+1, M, N)) + u_old(index('c', i+1, j+2, M, N)) + omega * sign(u_old(index('vx', i+1, j+2, M, N))) * (u_old(index('c', i+1, j+1, M, N)) - u_old(index('c', i+1, j+2, M, N))))/2;
            A(offset+(i-1)*M+j,index('vz', i+1, j+1, M, N)) = Ra*dt/hz(i+1) * (u_old(index('c', i+1, j+1, M, N)) + u_old(index('c', i, j+1, M, N)) + omega * sign(u_old(index('vz', i+1, j+1, M, N))) * (u_old(index('c', i+1, j+1, M, N)) - u_old(index('c', i, j+1, M, N))))/2;
            A(offset+(i-1)*M+j,index('vz', i+2, j+1, M, N)) = -Ra*dt/hz(i+1) * (u_old(index('c', i+2, j+1, M, N)) + u_old(index('c', i+1, j+2, M, N)) + omega * sign(u_old(index('vz', i+2, j+1, M, N))) * (u_old(index('c', i+2, j+1, M, N)) - u_old(index('c', i+1, j+2, M, N))))/2;   
        end
    end

    offset = offset + N*M;
    for j=1:M
        % Boundary conditions at the bottom
        A(offset+j,index('c', N+1, j+1, M, N)) = hz(N+2)/(hz(N+1)+hz(N+2));
        A(offset+j,index('c', N+2, j+1, M, N)) = hz(N+1)/(hz(N+1)+hz(N+2));
        A(offset+M+j, index('vz', N+2, j+1, M, N)) = 1;
                
        % Boundary conditions at the top
        A(offset+2*M+j, index('c', 1, j+1, M, N)) = hz(2)/(hz(1)+hz(2)) - 2/(hz(1)+hz(2));
        A(offset+2*M+j, index('c', 2, j+1, M, N)) = hz(1)/(hz(1)+hz(2)) + 2/(hz(1)+hz(2));
        A(offset+3*M+j, index('vz', 2, j+1, M, N)) = 1;
    end

    offset = offset + 4*M;
    for i=1:N
        % Periodic boundary conditions on the left and right side
        A(offset+i,index('c', i+1, 1, M, N)) = 1;
        A(offset+i,index('c', i+1, 2, M, N)) = 1;
        A(offset+i,index('c', i+1, M+1, M, N)) = -1;
        A(offset+i,index('c', i+1, M+2, M, N)) = -1;
        
        A(offset+N+i,index('c', i+1, 2, M, N)) = 1;
        A(offset+N+i,index('c', i+1, 1, M, N)) = -1;
        A(offset+N+i,index('c', i+1, M+2, M, N)) = -1;
        A(offset+N+i,index('c', i+1, M+1, M, N)) = 1;
        
        A(offset+2*N+i,index('p', i+1, 1, M, N)) = 1;
        A(offset+2*N+i,index('p', i+1, 2, M, N)) = 1;
        A(offset+2*N+i,index('p', i+1, M+1, M, N)) = -1;
        A(offset+2*N+i,index('p', i+1, M+2, M, N)) = -1;
        
        A(offset+3*N+i,index('p', i+1, 2, M, N)) = 1;
        A(offset+3*N+i,index('p', i+1, 1, M, N)) = -1;
        A(offset+3*N+i,index('p', i+1, M+2, M, N)) = -1;
        A(offset+3*N+i,index('p', i+1, M+1, M, N)) = 1;
        
        A(offset+4*N+i,index('vx', i+1, 2, M, N)) = 1;
        A(offset+4*N+i,index('vx', i+1, M+2, M, N)) = -1;
                
        if i>1
            A(offset+5*N+i-1,index('vz', i+1, 1, M, N)) = 1;
            A(offset+5*N+i-1,index('vz', i+1, 2, M, N)) = 1;
            A(offset+5*N+i-1,index('vz', i+1, M+2, M, N)) = -1;
            A(offset+5*N+i-1,index('vz', i+1, M+1, M, N)) = -1;
            
            A(offset+5*N+N-1+i-1,index('vz', i+1, 2, M, N)) = 1;
            A(offset+5*N+N-1+i-1,index('vz', i+1, 1, M, N)) = -1;
            A(offset+5*N+N-1+i-1,index('vz', i+1, M+2, M, N)) = -1;
            A(offset+5*N+N-1+i-1,index('vz', i+1, M+1, M, N)) = 1;
        end
    end
end

function u = initialize_state_cosine(M, N, alpha, a, hz, amplitude)
% Return a zero state vector, where the cosine perturbation seed has been added
% to the salt concentration
%
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - alpha: dimensionless height
% - a: perturbation wavenumber
% - hz: vector containing vertical cell sizes
% - amplitude: perturbation amplitude

    % Initialize the state vector to zero
    u = zeros(4*N*M+3*M+7*N-2,1);

    lambda = 2*pi/a;

    % x-coordinates of the cell middle points
    x = linspace(-lambda/(2*M),lambda+lambda/(2*M),M+2);

    % z-coordinates of the cell middle points
    z = zeros(N+2,1);
    for i=1:N+2
        z(i) = sum(hz(1:end-1)) - sum(hz(1:i-1)) - hz(i)/2;
    end

    inx = (z >= 0.75*alpha);

    % Add a cosine to the salinity in the top 25% of the domain
    for i = 1:nnz(inx)
        for j = 1:M+2
            if (i ~= 1 && i ~= N+2) || (j > 1 && j < M+2)
                u(index('c',i,j,M,N)) = amplitude * cos(a*x(j));
            end
        end
    end
end

function u = add_cosine(u, M, N, alpha, a, hz, amplitude)
% Add a cosine perturbation to the current salt concentration
%
% - u: current state vector
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - alpha: dimensionless height
% - a: perturbation wavenumber
% - hz: vector containing vertical cell sizes
% - amplitude: perturbation amplitude

    lambda = 2*pi/a;

    % x-coordinates of the cell middle points
    x = linspace(-lambda/(2*M),lambda+lambda/(2*M),M+2);

    % z-coordinates of the cell middle points
    z = zeros(N+2,1);
    for i=1:N+2
        z(i) = sum(hz(1:end-1)) - sum(hz(1:i-1)) - hz(i)/2;
    end

    inx = (z >= 0.75*alpha);

    % Add a cosine to the salinity in the top 25% of the domain
    for i = 1:nnz(inx)
        for j = 1:M+2
            if (i ~= 1 && i ~= N+2) || (j > 1 && j < M+2)
                u(index('c',i,j,M,N)) = u(index('c',i,j,M,N)) + amplitude * cos(a*x(j));
            end
        end
    end
end

function u = initialize_state_custom(M, N, alpha, Ra, a, tc, hz, amplitude)
% Return a zero state vector, where the customized perturbation seed (as predicted
% by the Galerkin projection) has been added to the salt concentration.
% Only works when a nonzero perturbation exists for the given parameters
% alpha, Ra, a and t_c
%
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - alpha: dimensionless height
% - Ra: Rayleigh number
% - a: perturbation wavenumber
% - tc: time of onset as predicted by the Cheybshev-Galerkin method for other given paramerts
% - hz: vector containing vertical cell sizes
% - amplitude: perturbation amplitude

    % Initialize the state vector to zero
    u = zeros(4*N*M+3*M+7*N-2,1);

    lambda = 2*pi/a;

    % x-coordinates of the cell middle points
    x = linspace(-lambda/(2*M),lambda+lambda/(2*M),M+2);

    % z-coordinates of the cell middle points
    z = zeros(N+2,1);
    for i=1:N+2
        z(i) = sum(hz(1:end-1)) - sum(hz(1:i-1)) - hz(i)/2;
    end

    % Get the salt perturbation profile
    perturb_c = Chebyshev_basis.salt_perturbation(alpha, a, Ra, tc);

    % Add the perturbation seed to the salinity
    for i = 1:N+2
        for j = 1:M+2
            if (i ~= 1 && i ~= N+2) || (j > 1 && j < M+2)
                u(index('c',i,j,M,N)) = u(index('c',i,j,M,N)) + amplitude * cos(a*x(j)) * perturb_c(z(i));
            end
        end
    end
end