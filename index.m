function ix = index(variable, i, j, M, N)
% Return the index of quantity 'c','p','vx' or 'vz' at position (i,j) in
% the state vector of the simulation
%
% - variable: 'c', 'p', 'vx' or 'vz'
% - i: vertical index
% - j: lateral index
% - M: amount of cells in lateral direction
% - N: amount of cells in vertical direction

   ix = 0;
   offset = 0;
   if variable == 'c'
       if i == 1 && j > 1 && j < M+2
           ix = j-1;
       elseif i > 1 && i < N+2 && j < M+3
           ix = M + (i-1)*(M+2) - (M+2-j);
       elseif i == N+2 && j > 1 && j < M+2
           ix = (M+2)*(N+2)-4 - (M+1-j);
       end          
   elseif variable == 'p'
       offset = (M+2)*(N+2)-4;
       if i > 1 && i < N+2 && j < M+3
           ix = (i-1)*(M+2) - (M+2-j);
       end
   elseif variable == 'vx'
       offset = (M+2)*(N+2)-4+(M+2)*N;
       if i > 1 && i < N+2 && j > 1 && j < M+3
           ix = (i-1) * (M+1) - (M+2-j);
       end
   elseif variable == 'vz'
       offset = (M+2)*(N+2)-4+N*(M+2)+N*(M+1);
       if i == 2 && j > 1 && j < M+2
           ix = j-1;
       elseif i > 2 && i < N+2 && j < M+3
           ix = M + (i-2)*(M+2) - (M+2-j);
       elseif i == N+2 && j > 1 && j < M+2
           ix = (M+2)*(N+1)-4 -(M+1-j);
       end
   end
   if ix ~= 0
       ix = ix + offset;
   end
end