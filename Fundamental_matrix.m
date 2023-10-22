classdef Fundamental_matrix
    methods (Static)

        function y = determinant_function(alpha, a, t, c_S)
        % Returns the determinant of the matrix M as defined in the fundamental matrix method 
        % as a function of the Rayleigh number for the bisection
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - t: time of the ground state
        % - c_S: Ground_state object for given alpha

            y = @(Ra) Fundamental_matrix.get_determinant(alpha, a, Ra, t, c_S);
        end
        
        function d = get_determinant(alpha, a, Ra, t, c_S)
        % Returns the determinant of the matrix M as defined in the fundamental matrix method 
        % for a given set of parameters
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - Ra: Rayleigh number
        % - t: time of the ground state
        % - c_S: Ground_state object for given alpha


            % Construct the matrices describing the boundary conditions of
            % the perturbation problem
            R1 = zeros(4,4);
            R1(1,1) = 1; 
            R1(3,3) = 1;
            R2 = zeros(4,4);
            R2(2,1) = 1;
            R2(4,3) = 1;
            R2(4,4) = -1;

            % Get the fundamental matrix \Phi(\alpha;0)
            Phi = Fundamental_matrix.get_fundamental_matrix(alpha, a, Ra, t, c_S);

            % Calculate the determinant of M
            d = det(R1 + R2*Phi);
        end
        
        function Phi = get_fundamental_matrix(alpha, a, Ra, t, c_S)
        % Returns the fundamental matrix Phi(alpha;0) of the perturbation system 
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - Ra: Rayleigh number
        % - t: time of the ground state
        % - c_S: Ground_state object for given alpha


            N = size(c_S.spatial_derivatives_at_zero,2);

            % Initialize cell arrays to save the power series coefficient 
            % matrices of A(z) and Phi(z;0)
            A_k = cell(N,1);
            A_k{1} = [0,1,0,0; a^2,0,a^2,0; 0,0,0,1; Ra*c_S.get_derivative_at_zero(t,1),0,a^2,1];
            Phi_k = cell(N,1);
            Phi_k{1} = eye(4);
            Phi = Phi_k{1};
            
            % Calculate the coefficient matrices A_k of the system matrix A(z)
            for k = 2:N
                T = zeros(4,4);
                T(4,1) = Ra * c_S.get_derivative_at_zero(t, k) / factorial(k-1);
                A_k{k} = T;
            end
            
            % Calculate the coefficient matrices Phi_k of the fundamental
            % matrix Phi(z;0)
            for k = 1:N
                T = zeros(4,4);
                for j=1:k
                    T = T + A_k{k-j+1} * Phi_k{j};
                end
                Phi_k{k+1} = (1/k) * T;

                % Evaluate the power series of the fundamental matrix at alpha
                Phi = Phi + alpha^k * Phi_k{k+1};
            end            
        end
        
        function A = get_system_matrix(alpha, a, Ra, t, c_S, z)
        % Returns the system matrix A(z) of the perturbation system for the
        % given parameter combination
        %
        % - alpha: dimensionless height
        % - a: perturbation wavenumber
        % - Ra: Rayleigh number
        % - t: time of the ground state
        % - c_S: Ground_state object for given alpha
        % - z: height

            A = [0, 1, 0, 0; a^2, 0, a^2, 0; 0, 0, 0, 1; Ra*c_S.get_derivative(z,t), 0, a^2, 1];
        end
    end
end