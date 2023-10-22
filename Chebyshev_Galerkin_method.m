% Script that calculates the eigenvalue Ra(alpha, t, a) of the perturbation
% problem for a given set of parameter combinations with the Chebyshev-Galerkin 
% method. The results are saved in the 4 x n_alpha x n_time x n_a array 
% Ra_alpha_t_a, which contains Ra in the first entry, alpha in the second entry, 
% t in the third entry and a in the fourth entry. 
% In this script, only the parameter grids should be adjusted.

% Parameters grids for which Ra(alpha, t, a) should be calculated via
% the Chebyshev-Galerkin approach.
% Enter one or more values of interest for every variable (Ra is computed 
% for every possible combination of the values)
time = [0.2];
alpha = [1];
a = [0.01];

% Number of Chebyshev polynomials to be used for omega
m = 30;
% Number of Chebyshev polynomials to be used for chi
n = 30;

% Sizes of the grids
n_time = max(size(time));
n_alpha = max(size(alpha));
n_a = max(size(a));

% Array to save the Rayleigh number as function of alpha, time and a
Ra_alpha_t_a = zeros(4,n_alpha,n_time,n_a);

% Iterate over all alphas
for i = 1:n_alpha
    
    % Calculate the ground-state salinity
    ground_state = Ground_state(alpha(i), 100, @(z) 0);

    % Iterate over all times
    for j = 1:n_time
        
        % Get the first spatial derivative of the ground-state salinity
        ground_state_der_val = ground_state.get_derivative(linspace(0,alpha(i),1000), time(j));
        ground_state_der = @(z) interp1(linspace(0,alpha(i),1000),ground_state_der_val, z);
        
        % Iterate over all wavenumbers
        for k = 1:n_a

            % Get the matrix M describing the linear equation system
            M = Chebyshev_basis.get_matrix_M(alpha(i), a(k), 1, m, n, ground_state_der);
            U = M(1:m,1:m);
            V = M(1:m,m+1:end);
            W = M(m+1:end,1:m);
            X = M(m+1:end,m+1:end);

            % Solve the generalized eigenvalue problem
            H = linsolve(X,W);
            eigs = eig(U, V*H);
                
            % Filter unphysical eigenvalues
            index = (eigs >= 0) & (~isinf(eigs)) & (imag(eigs) == 0);
            
            % If there is a physical eigenvalue, set the Rayleigh number to
            % the smallest, else set it to nan
            if nnz(index) > 0
                Ra_alpha_t_a(1,i,j,k) = min(eigs(index));
            else
                Ra_alpha_t_a(1,i,j,k) = nan;
            end

            Ra_alpha_t_a(2,i,j,k) = alpha(i);
            Ra_alpha_t_a(3,i,j,k) = time(j);
            Ra_alpha_t_a(4,i,j,k) = a(k);
        end
    end
    Ra_alpha_t_a(:,i,:,:)
end