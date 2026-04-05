%% Q2 Take Home Final 1D
gam = 3e-5;
Nx = 128;
Lx = 1;
TT = 0.05;
dt = 1e-6;

Nt = fix(TT/dt);
dx = Lx/Nx;
x = (0:Nx-1)'*dx;

Dxx = toeplitz([-2 1 zeros(1,Nx-3) 1])/dx/dx;
Dxx = sparse(Dxx);

c = rand(Nx,1) > 1/2;

cs = zeros(Nx, Nt); % save every time step

for nt = 1:Nt
    % calculate mu from equations [2] and [3]
    W = (c - 1) .* c .* (c - 0.5);   % W(c) from equation [3]
    mu = W - gam * (Dxx * c);          % mu from equation [2]
    
    dc = Dxx * mu;                      % from equation [1]
    c = c + dt * dc;                    % forward Euler update rule
    
    % plot feedback
    if(rem(nt, fix(Nt/100)) == 0)
        plot(x, c)
        drawnow;
        disp([nt/100 mean(c(:))]);
    end
    
    cs(:, nt) = c;                      % save results
end

% Final plot
imagesc([0 Lx],[0 TT], cs');
xlabel('Space');
ylabel('Time');
colormap(jet(256));