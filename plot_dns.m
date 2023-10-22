function plot_dns(u, alpha, a, M, N)
% Plots the result of the direct numerical simulation
%
% - u: state vector of the simulation result to be plotted
% - alpha: dimensionless height
% - a: perturbation wavenumber used in the simulation
% - M: amount of lateral cells used in the simulation
% - N: amount of vertical cells used in the simulation

    lambda = 2*pi/a;
    [X,Z] = meshgrid(linspace(-lambda/(2*M),lambda*(2*M+1)/(2*M),M+2),linspace(alpha*(2*N+1)/(2*N),-alpha/(2*N),N+2));
    [X__,Z__] = meshgrid(linspace(0,2*lambda,99),linspace(0,alpha,50));
    [X_,Z_] = meshgrid(linspace(0,lambda,50),linspace(0,alpha,50));
    [X_vx,Z_vx] = meshgrid(linspace(0,lambda,M+1),linspace(-alpha/(2*N),alpha*(2*N+1)/(2*N),N+2));
    [X_vz,Z_vz] = meshgrid(linspace(-lambda/(2*M),lambda*(2*M+1)/(2*M),M+2),linspace(0,alpha,N+1));

    u1 = transform(M, N, u);

    u1.c(1,1) = u1.c(1,2);
    u1.c(1,end) = u1.c(1,end-1);
    u1.c(end,1) = u1.c(end,2);
    u1.c(end,end) = u1.c(end,end-1);
    u1.c = flip(u1.c,2);
    c_ = @(x,z) interp2(X,Z,u1.c,x,z);
    c__ = c_(X_,Z_);

    figure;
    pcolor(X__,Z__,[c__,c__(:,2:end)])

    u1.vx(1,:) = u1.vx(2,:);
    u1.vx(end,:) = u1.vx(end-1,:);
    u1.vx = u1.vx(:,2:end-1);
    u1.vx = flip(flip(u1.vx,2),1);
    u1.vx = -u1.vx;
    vx_ = @(x,z) interp2(X_vx,Z_vx,u1.vx,x,z);
    vx__ = vx_(X_,Z_);

    u1.vz(2,1)=u1.vz(2,2);
    u1.vz(end-1,1)=u1.vz(end-1,2);
    u1.vz(2,end)=u1.vz(2,end-1);
    u1.vz(end-1,end)=u1.vz(end-1,end-1);
    u1.vz = u1.vz(2:end-1,:);
    u1.vz = flip(flip(u1.vz,2),1);
    vz_ = @(x,z) interp2(X_vz,Z_vz,u1.vz,x,z);
    vz__ = vz_(X_,Z_);
    
    shading interp
    colormap('turbo')
    hold on 
    brighten(0.1)
    q = streamslice(X__,Z__,[vx__,vx__(:,2:end)],[vz__,vz__(:,2:end)],'cubic');
    set(q,'color','k')
    set(q,'Linewidth',1.2)
    cbar = colorbar;
    set(gca, 'fontname', 'times new roman')
    set(gca, 'fontsize', 15)
    xlabel('$\hat{x}$','interpreter','latex','fontsize',20)
    ylabel('$\hat{z}$','interpreter','latex','fontsize',20)
    cbar.Label.Interpreter = 'latex';
    cbar.Label.String = '$\hat{c}$';
    cbar.Label.FontSize = 20;
    Ra = 1/u1.vz(1,2);
    title(sprintf('$\\alpha = %0.1f, \\mathrm{Ra} = %0.1f, \\hat{a} = %0.2f$',alpha,Ra,a),'interpreter','latex','fontsize',22,'fontname','times new roman')
end