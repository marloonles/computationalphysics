%% 2D Cahn-Hilliard Simulator
gam = 3e-4;
Nx = 128;
Ny = 128;
Lx = 1;
Ly = 1;
TT = 0.005;
dt = 1e-7;

%% Calculated parameters
Nt = fix(TT/dt);
dx = Lx/Nx;
dy = Ly/Ny;

Dxx = toeplitz([-2 1 zeros(1,Nx-3) 1])/dx/dx;
Dxx = sparse(Dxx);

Dyy = toeplitz([-2 1 zeros(1,Ny-3) 1])/dy/dy;
Dyy = sparse(Dyy);

%% Initial conditions
c = rand(Nx, Ny) > 1/2;

%% Main loop
for nt = 1:Nt
    W = (c - 1) .* c .* (c - 0.5);
    mu = W - gam * (Dxx*c + c*Dyy);
    dc = Dxx*mu + mu*Dyy;
    c = c + dt*dc;

    if(rem(nt, fix(Nt/200)) == 0)
        imagesc([0 Ly],[0 Lx], c);
        axis('image');
        colormap(jet(256));
        drawnow;
        disp([nt/100 mean(c(:))]);
    end
end