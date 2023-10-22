% Script that calculates the eigenvalue Ra(alpha, t, a) of the perturbation
% problem for a given set of parameter combinations with the fundamental matrix
% method. The results are saved in the 4 x n_alpha x n_time x n_a array 
% Ra_alpha_t_a, which contains Ra in the first entry, alpha in the second 
% entry, t in the third entry and a in the fourth entry.
% In this script, only the parameter grids and the grid of Rayleigh numbers 
% used to find the eigenvalue should be adjusted.

% Parameters grids for which Ra(alpha, t, a) should be calculated via
% the fundamental matrix method.
% Enter one or more values of interest for every variable (Ra is computed
% for every possible combination of the values)
time = [0.2];
alpha = [1];
a = [0.01];

% Grid of Rayleigh numbers to search a change of sign 
% of the determinant of M to start a bisection.
% Ra(alpha, t, a) has to be between the smallest and largest value in this
% array, otherwise the bisection won't find an eigenvalue!
% (This array is used to find a change of sign of the determinant of M from
% equation (45))
Ra = logspace(-8, 6, 40);

% Sizes of the grids
n_time = max(size(time));
n_alpha = max(size(alpha));
n_a = max(size(a));
n_r = max(size(Ra));

% Array to save the Rayleigh number as function of alpha, time and a
Ra_alpha_t_a = zeros(4,n_alpha,n_time,n_a);

% Iterate over all alphas
for i = 1:n_alpha

    % Calculate the ground-state salinity
    c_S = Ground_state(alpha(i), 100, @(z) 0);
    c_S.get_spatial_derivatives_at_zero(100);

    % Iterate over all times
    for j = 1:n_time

        % Iterate over all wavenumbers
        for k = 1:n_a
            
            % Find Rayleigh numbers such that the determinant of M changes
            % sign in order to start the bisection
            interval = zeros(2,1);
            found = false;
            l = 1;
            dstart = Fundamental_matrix.get_determinant(alpha(i), a(k), Ra(1), time(j), c_S);
            interval(1) = Ra(1);
            dend = 0;
            
            while l < n_r && ~found

                l = l+1;
                dend = Fundamental_matrix.get_determinant(alpha(i), a(k), Ra(l), time(j), c_S);
                
                if sign(dstart) == sign(dend) | isinf(dstart) | isinf(dend) | isnan(dstart) | isnan(dend)
                    dstart = dend;
                    interval(1) = Ra(l);
                else
                    found = true;
                    interval(2) = Ra(l);
                end
            end
                        
            % Employ bisection if change of sign of determinant was found
            if found
               Ra_alpha_t_a(1,i,j,k) = fzero(Fundamental_matrix.determinant_function(alpha(i), a(k), time(j), c_S), interval);  
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