classdef Ground_state < handle
% Ground_state Class which handles the calculation of the ground-state salt concentration

    properties
        spatial_functions           % Sturm-Liouville eigenfunctions
        temporal_functions          % Temporal eigenfunctions
        initial_condition           % Function handle describing the initial condition for the salinity
        initial_coefficients        % Initial coefficients a_n for every eigenfunction
        spatial_derivatives         % First derivatives of the Sturm-Liouville eigenfunctions
        spatial_derivatives_at_zero % Higher-order derivatives of the Sturm-Liouville eigenfunctions evaluated at zero
        alpha                       % Dimensionless height
        delta                       % Wavenumbers appearing in the Sturm-Liouville eigenfunctions
        eps                         % Epsilon appearing in the first Sturm-Liouville eigenfunction
        n                           % Truncation of the Sturm-Liouville Fourier series
    end
    
    methods
        function obj = Ground_state(alpha, N, phi)
        % Constructor of the Ground_state
        %
        % - alpha: dimensionless height
        % - N: amount of Sturm-Liouville eigenfunctions to be used
        % - phi: Function handle of the initial condition (usually: @(z) 0)

            obj.alpha = alpha;
            obj.n = N;

            if obj.alpha > 2
                obj.eps = obj.get_epsilon();
            end

            obj.delta = obj.get_delta();
            obj.initial_condition = @(z) phi(z) - exp(z) + 1;
            obj.spatial_functions = obj.get_spatial_functions();
            obj.spatial_derivatives = obj.get_spatial_derivatives();
            obj.temporal_functions = obj.get_temporal_functions();
            obj.initial_coefficients = obj.get_initial_coefficients();
        end
        
        function c = get_solution(obj, z, t)
        % Return the value of the ground-state salt concentration at
        % heights z and time t
        %
        % - z: vector of z values
        % - t: time at which ground state is evaluated

            q = zeros(size(z));
            for i=1:obj.n
                q = q + obj.initial_coefficients(i) .* obj.temporal_functions{i}(t) .* obj.spatial_functions{i}(z);
            end
            c = q + exp(z) -1;
        end
        
        function c = get_derivative(obj, z, t)
        % Return the value of the first spatial derivative of the ground-state 
        % salt concentration at points z and time t
        % 
        % - z: vector of z values
        % - t: time at which derivative is evaluated
         
            q = zeros(size(z));
            for i=1:obj.n
                q = q + obj.initial_coefficients(i) .* obj.temporal_functions{i}(t) .* obj.spatial_derivatives{i}(z);
            end
            c = q + exp(z);
        end
        
        function delta = get_delta(obj)
        % Return the wavenumbers delta of the Sturm-Liouville
        % eigenfunctions

            delta = zeros(obj.n,1);
            funn = @(x) 0.5*tan(x*obj.alpha)-x;
            if obj.alpha < 2
                x1 = 10^(-10);
                x2 = pi/(2*obj.alpha);
                delta(1) = fzero(funn,[x1,x2]);
                for i = 1:obj.n-1
                    x1 = (2*i-1)*pi/(2*obj.alpha) + 10^(-10);
                    x2 = (2*i+1)*pi/(2*obj.alpha) - 10^(-10);
                    delta(i+1) = fzero(funn,[x1,x2]);
                end
            else
                for i = 1:obj.n
                    x1 = (2*i-1)*pi/(2*obj.alpha) + 10^(-10);
                    x2 = (2*i+1)*pi/(2*obj.alpha) - 10^(-10);
                    delta(i) = fzero(funn,[x1,x2]);
                end
            end
        end
        
        function res = get_epsilon(obj)
        % Return the epsilon parameter appearing in the first
        % Sturm-Liouville eigenfunction for alpha > 2

            fx = @(x) exp(2*x*obj.alpha) * ((0.5-x)/(0.5+x)) -1;
            res = fzero(fx,[10^-10,0.5 + 10^-10]);
        end
        
        function spatial_functions = get_spatial_functions(obj)
        % Return a cell array of function handles of the first obj.n
        % Sturm-Liouville eigenfunctions

            spatial_functions = cell(obj.n,1);
            first = 0;
            if obj.alpha == 2
                spatial_functions{1} = @(z) z.*exp(0.5*z);
                first = 1;
            end
            if obj.alpha > 2
                r1 = 0.5-obj.eps;
                r2 = 0.5+obj.eps;
                spatial_functions{1} = @(z) -exp(r1*z+2*obj.eps*obj.alpha) *(r1/r2) + exp(r2*z);
                first = 1;
            end
            if first
                for k=2:obj.n
                    spatial_functions{k} = @(z) exp((1/2)*z) .* sin(obj.delta(k-1)*z);
                end
            else
                for k=1:obj.n
                    spatial_functions{k} = @(z) exp((1/2)*z) .* sin(obj.delta(k)*z);
                end
            end
        end
        
        function spatial_derivatives = get_spatial_derivatives(obj)
        % Return a cell array of function handles of the first 
        % spatial derivative of the first obj.n Sturm-Liouville 
        % eigenfunctions

            spatial_derivatives = cell(obj.n,1);
            first = 0;
            if obj.alpha == 2
                spatial_derivatives{1} = @(z) (0.5.*z + 1) .* exp(0.5*z);
                first = 1;
            end
            if obj.alpha > 2
                r1 = 0.5-obj.eps;
                r2 = 0.5+obj.eps;
                spatial_derivatives{1} = @(z) -r1*exp(r1*z+2*obj.eps*obj.alpha) * (r1/r2) + r2*exp(r2*z);
                first = 1;
            end
            if first
                for k=2:obj.n
                    spatial_derivatives{k} = @(z) exp((1/2)*z) .* (0.5*sin(obj.delta(k-1)*z) + obj.delta(k-1)*cos(obj.delta(k-1).*z));
                end
            else
                for k=1:obj.n
                    spatial_derivatives{k} = @(z) exp((1/2)*z) .* (0.5*sin(obj.delta(k)*z) + obj.delta(k)*cos(obj.delta(k).*z));
                end
            end
        end
        
        function temporal_functions = get_temporal_functions(obj)
        % Return a cell array of function handles of the first obj.n
        % temporal eigenfunctions Gamma_n

            temporal_functions = cell(obj.n,1);
            if obj.alpha >= 2
                if obj.alpha == 2
                    temporal_functions{1} = @(t) exp(-0.25*t);
                end
                if obj.alpha > 2
                    temporal_functions{1} = @(t) exp(-(0.25-obj.eps^2)*t);
                end
                for k=2:obj.n
                    temporal_functions{k} = @(t) exp(-(obj.delta(k-1)^2 + 0.25)*t);
                end
            else
                for k=1:obj.n
                    temporal_functions{k} = @(t) exp(-(obj.delta(k)^2 + 0.25)*t);
                end
            end
        end
         
        function an = get_initial_coefficients(obj)
        % Return the initial coefficients a_n of the corresponding
        % Sturm-Liouville eigenfunction X_n

           an = ones(obj.n,1);
           norms = obj.get_norms();
           for k=1:obj.n
               X = obj.spatial_functions{k};
               func = @(z) obj.initial_condition(z) .* X(z) .* exp(-z);
               an(k) = (1/norms(k)) * integral(func, 0, obj.alpha);
           end
        end
        
        function xn = get_norms(obj)
        % Return the norms of the Sturm-Liouville eigenfunctions with 
        % respect to the weighted inner product of L^2([0,alpha],exp(-z))

            xn = zeros(obj.n,1);
            if obj.alpha >= 2
                fz = @(z) obj.spatial_functions{1}(z).^2 .* exp(-z);
                xn(1) = integral(fz, 0, obj.alpha);
                for k=2:obj.n
                    xn(k) = 0.5*obj.alpha - sin(2*obj.alpha*obj.delta(k-1))/(4*obj.delta(k-1));
                end
            else
                for k=1:obj.n
                    xn(k) = 0.5*obj.alpha - sin(2*obj.alpha*obj.delta(k))/(4*obj.delta(k));
                end
            end
        end

        function der = get_derivative_at_zero(obj, t, k)
        % Return the value of the k-th spatial derivative of the
        % ground-state salt concentration at time t and z = 0
        %
        % - t: time at which the derivative is evaluated
        % - k: order of the derivative

            der = 1;
            for i = 1:obj.n
                der = der + obj.initial_coefficients(i) * obj.temporal_functions{i}(t) * obj.spatial_derivatives_at_zero(i,k);
            end
        end
        
        function der = get_spatial_derivatives_at_zero(obj, kmax)
        % Return an array containing the values of the first k spatial
        % derivative of the obj.n Sturm-Liouville eigenfunctions at z = 0
        %
        % - kmax: order of the highest derivative which is computed

            der = zeros(obj.n,kmax);
            
            if obj.alpha >= 2
                if obj.alpha == 2
                    der(1,1) = 1;
                    for k = 2:kmax
                        der(1,k) = 0.5*der(1,k-1) + (0.5)^(k-1); 
                    end
                elseif obj.alpha > 2
                    r1 = 0.5 + obj.eps;
                    r2 = 0.5 - obj.eps;
                    for k = 1:kmax
                        der(1,k) = r1^k - r2^k * exp(2*obj.eps*obj.alpha) * (r2/r1);
                    end
                end

                for i = 2:obj.n
                    for k = 1:kmax
                        der(i,k) = ((0.5 + obj.delta(i-1)*1i)^k - (0.5 - obj.delta(i-1)*1i)^k) / 2i;
                    end
                end

            else
                for i = 1:obj.n
                    for k = 1:kmax
                        der(i,k) = ((0.5 + obj.delta(i)*1i)^k - (0.5 - obj.delta(i)*1i)^k) / 2i;
                    end
                end
            end

            obj.spatial_derivatives_at_zero = der;
        end
    end
end