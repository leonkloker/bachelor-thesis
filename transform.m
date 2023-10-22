function u = transform(M, N, x)
% Transforms the state vector x of the simulation into a struct with 2-dimensional arrays for
% pressure, vertical- and lateral velocity and salinity
%
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction
% - x: state vector

    u = struct;
    u.c = zeros(N+2,M+2);
    u.p = zeros(N+2,M+2);
    u.vx = zeros(N+2,M+3);
    u.vz = zeros(N+3,M+2);
    
    u.vx(1,:) = NaN;
    u.vx(end,:) = NaN;
    u.vx(:,1) = NaN;
    u.vx(:,end) = NaN;
    u.vz(1,:) = NaN;
    u.vz(2,end) = NaN;
    u.vz(2,1) = NaN;
    u.vz(end,:) = NaN;
    u.vz(end-1,end) = NaN;
    u.vz(end-1,1) = NaN;
    u.c(1,1) = NaN;
    u.c(1,end) = NaN;
    u.c(end,1) = NaN;
    u.c(end,end) = NaN;
    u.p(end,:) = NaN;
    u.p(1,:) = NaN;
    
    indc = ~isnan(u.c);
    indp = ~isnan(u.p);
    indvx = ~isnan(u.vx);
    indvz = ~isnan(u.vz);
    
    ch = NaN((M+2)*(N+2),1);
    ph = NaN((M+2)*(N+2),1);
    vxh = NaN((M+2)*(N+2),1);
    vzh = NaN((M+2)*(N+2),1);
    
    for j=1:M+2
        for i=1:N+2
            if index('c', i, j, M, N) ~= 0
                ch((j-1)*(N+2)+i) = x(index('c', i, j, M, N));
            end
            if index('p', i, j, M, N) ~= 0
                ph((j-1)*(N+2)+i) = x(index('p', i, j, M, N));
            end
            if index('vx', i, j, M, N) ~= 0
                vxh((j-1)*(N+2)+i) = x(index('vx', i, j, M, N));
            end
            if index('vz', i, j, M, N) ~= 0
                vzh((j-1)*(N+2)+i) = x(index('vz', i, j, M, N));
            end
        end
    end
 
    u.c(indc) = ch(~isnan(ch));
    u.p(indp) = ph(~isnan(ph));
    u.vx(indvx) = vxh(~isnan(vxh));
    u.vz(indvz) = vzh(~isnan(vzh));
end