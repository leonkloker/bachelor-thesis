classdef Chebyshev_basis
    methods (Static)
        
        function M = get_matrix_M(alpha, a, Ra, m, n, dc_S)
        % Construct the matrix M describing the linear equation system
        % arising from the Galerkin projection of the perturbation problem
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - Ra: Rayleigh number
        % - m: number of Chebyshev polynomials to be used for velocity perturbation omega
        % - n: number of Chebyshev polynomials to be used for salt perturbation chi
        % - dc_S: function handle of the first spatial derivative of the ground state
            
            % Load all parameter-independent integral values from file
            values = load('Chebyshev_integrals.mat');
            
            % Construct first constituent matrix of M
            U = zeros(m,m);
            for i=1:m-2
                for j=1:m
                    U(i,j) = (2/alpha)*values.int_ddTT(i,j) - (a^2 *alpha/2)*values.int_TT(i,j);
                end
            end

            % Construct second constituent matrix of M
            V = zeros(m,n);
            for i=1:m-2
                for j=1:n
                    V(i,j) = values.int_TT(i,j);
                end
            end
            V = -(a^2 * alpha/2) * V;
            
            % Construct third constituent matrix of M
            W = zeros(n,m);
            for i=1:n-2
                for j=1:m
                    % This integral is calculated via Gaussian quadrature,
                    % since it can't be precomputed as it depends on dc_S
                    W(i,j) = Chebyshev_basis.gaussian_quadrature_func_p_p(alpha, i-1, j-1, @(z) dc_S(z));
                end
            end
            W = -Ra * (alpha/2) .* W;
            
            % Construct fourth constituent matrix of M
            X = zeros(n,n);
            for i=1:n-2
                for j=1:n
                    X(i,j) = (2/alpha)*values.int_ddTT(i,j) - values.int_dTT(i,j) - (a^2 * alpha/2)*values.int_TT(i,j);
                end
            end
            
            % Include the flow boundary conditions as additional equations
            % in U and V
            for j=1:m
                U(m-1,j) = (-1)^(j-1);
                U(m,j) = 1;
            end
           
            % Include the salt boundary conditions as additional equations     
            % in W and X
            for j=1:n
                X(n-1,j) = (-1)^(j-1);
                X(n,j) = 1 - (2/alpha)*(j-1)^2;
            end
            
            % Assemble M
            M = zeros(m+n,m+n);
            M(1:m,1:m) = U;
            M(1:m,m+1:end) = V;
            M(m+1:end,1:m) = W;
            M(m+1:end,m+1:end) = X;            
        end
        
        function chi = salt_perturbation(alpha, a, Ra, t)
        % Calculate the profile of the vertical salt perturbation and
        % return it as a function handle
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - Ra: Rayleigh number
        % - t: time of the system, such that a nonzero perturbation exists
            
            % Default amount of Chebyshev polynomials used for the Galerkin
            % projection
            m = 30;
            n = 30;
            
            % Calculate ground-state salinity and its derivative
            ground_state = Ground_state(alpha, 100, @(z) 0);
            ground_state_der_val = ground_state.get_derivative(linspace(0,alpha,1000), t);
            ground_state_der = @(z) interp1(linspace(0,alpha,1000),ground_state_der_val, z);
            
            % Get matrix M of the Galerkin projection
            M = Chebyshev_basis.get_matrix_M(alpha, a, Ra, m, n, @(z) ground_state_der(z));
            
            % Calculate the nullspace of M to get admissible perturbations
            coeff = null(M);
            
            % If the nullspace is empty, the matrix was non-degenerate and
            % no neutral perturbation exists
            if min(size(coeff)) == 0 || sum(isnan(coeff)) ~= 0
                error('There is no non-zero perturbation for these parameters!');
                return;
            end

            z = linspace(0, alpha, 1000);
            
            % Calculate the shape of the predicted salt perturbation
            chi = zeros(size(z));
            for s = 1:n
                chi = chi + coeff(m+s)*chebyshevT(s-1,(2/alpha)*z-1);
            end
               
            % Normalize the salt perturbation and return it as a function
            % handle
            chi = @(x) interp1(z, chi./max(abs(chi)), x, 'linear', 'extrap');
        end
        
        function integral = gaussian_quadrature_func_p_p(alpha, i, j, f)
        % Calculate the integral of T_i(t) * T_j(t) * f((alpha/2)*(t+1)) from -1 to 1 via Gaussian quadrature 
        %
        % - alpha: dimensionless height
        % - i: first Chebyshev polynomial to be used
        % - j: second Chebyshev polynomial to be used
        % - f: function handle of the function to be integrated
  
           % Load the evaluation points, quadrature weights and Cheybshev
           % polynomial values at the evaluation points from a file
           quadrature = load('Gaussian_quadrature.mat');

           weights = quadrature.weights;
           values = quadrature.chebyshev_values;
           points = quadrature.points;
            
           % Calculate the integral with the given quadrature weights and
           % points
           integral = 0;
           for k = 1:max(size(points))
              integral = integral + weights(k) * values(k,j+1) * values(k,i+1) * f((alpha/2)*(points(k)+1));
           end
        end
        
        function precalculate_chebyshev_integrals()
        % Calculate the integrals of the Chebyshev polynomials with each
        % other, that are needed for the Galerkin projection and save their
        % values in a file

            % Maximal amount of Chebyshev polynomials to be used for the Galerkin
            % projection
            n = 30;
            
            % Arrays which save the values of the calculated integrals
            int_TT = zeros(n);
            int_dTT = zeros(n);
            int_ddTT = zeros(n);

            syms z;

            for i=1:n
                for j=1:n
                    T_j = sym2poly(chebyshevT(j-1,z));
                    T_i = sym2poly(chebyshevT(i-1,z));
                    p = conv(T_i,T_j);
        
                    % int_-1^1 T_j-1 * T_i-1 dx
                    int_TT(i,j) = diff(polyval(polyint(p),[-1 1]));
                    
                    dT_j = polyder(T_j);
                    p = conv(T_i,dT_j);
        
                    % int_-1^1 dT_j-1 * T_i-1 dx 
                    int_dTT(i,j) = diff(polyval(polyint(p),[-1 1]));
        
                    ddT_j = polyder(dT_j);
                    p = conv(T_i,ddT_j);
        
                    % int_-1^1 d^2T_j-1 * T_i-1 dx
                    int_ddTT(i,j) = diff(polyval(polyint(p),[-1 1]));
                end
            end

            % Save the integrals
            save('Chebyshev_integrals.mat','int_TT');
            save('Chebyshev_integrals.mat','int_dTT','-append');
            save('Chebyshev_integrals.mat','int_ddTT','-append');
        end

        function precalculate_gaussian_quadrature()
        % Evaluate the Chebyshev polynomials at the points of the Gaussian
        % quadrature on [-1,1] and save them in a file together with quadrature points
        % and weights
            
            % Maximal amount of Chebyshev polynomials to be used for the Galerkin
            % projection
            n = 30;

            % Gaussian quadrature with 40 collocation points and corresponding weights
            chebyshev_values = zeros(40,n);

            % Get the points and weights of the Gaussian quadrature
            [points,weights] = lgwt(40,-1,1);

            for i = 1:n
                for j = 1:40
                    
                    % Evaluate the Chebyshev i-1-th Chebyshev polynomial at the j-th
                    % collocation point
                    chebyshev_values(j,i) = chebyshevT(i-1,points(j));  
                end
            end

            % Save the Gaussian quadrature data
            save('Gaussian_quadrature.mat','chebyshev_values');
            save('Gaussian_quadrature.mat','weights','-append');
            save('Gaussian_quadrature.mat','points','-append');
        end
    end
end