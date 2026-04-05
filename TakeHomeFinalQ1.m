N=5; % Number of nodes
a=3/8; % Location of forcing delta function

dx=1/(N-1); % grid spacing Delta

Kloc = (1/30)*[ 4    3   -1   -3;
                3   36    3  -36;
               -1    3    4   -3;
               -3  -36   -3   36];

Kloc = Kloc/dx;

%% Assemble global K (10x10)
K=zeros(2*N);
for i=0:N-2
    idx=(2*i+1):(2*i+4);
    K(idx,idx)=K(idx,idx)+Kloc;
end

%% Build f (10x1)
S=@(nu) nu.*(abs(nu)-1).^2.*(abs(nu)<1);
H=@(nu) (2*abs(nu)+1).*(abs(nu)-1).^2.*(abs(nu)<1);

f=zeros(2*N,1);
for n=0:N-1
    nu=a/dx-n;
    f(2*n+1)=S(nu)*dx;
    f(2*n+2)=H(nu)*dx;
end

%% Apply boundary conditions
% u'(0)=0 => u^S_0=0:
% u(1)=0  => u^H_4=0:
dofs=2:2*N-1;
Kbc=K(dofs,dofs);
fbc=f(dofs);

%% Solve for u
u_inner=Kbc\fbc;

% reinsert boundary values
u=[0; u_inner; 0]; % u^S_0=0 prepended, u^H_4=0 appended

%% Finite differences (free-fixed)
e = ones(N-1,1);
Dxx = -(diag(e(1:end-1),1) - 2*diag(e) + diag(e(1:end-1),-1))/dx^2;
Dxx(1,1) = 1/dx^2;  % Neumann BC ghost point

fFD = zeros(N-1,1);
fFD(2) = 1/2;
fFD(3) = 1/2;

u2 = Dxx\fFD;
u2 = [u2; 0];  % append u(1)=0

%% Evaluate FEM solution
x=0:dx/50:1; % locations to find U(x)

U = evalFEM(x,u,dx) % calculate U(x) from x, u, and dx

ue = (1-x)/4.*(x>a) + (1-a)/4.*(x<=a) % exact solution

%% Plot results
h=plot(x,U,(0:N-1)*dx,u2,'ro--',x,ue,'k--');

%% Make it pretty
set(h,'linewidth',3,'markersize',15);
set(gca,'fontsize',20);
xlabel('$x$','interp','latex');
ylabel('$u(x)$','interp','latex');
legend('FEM','FD','Exact');