clear variables
L = 1;
Nx = 100;
gridx = linspace(0,L,Nx);
Nv = 10;
Nt = 50;
bins = 20;
T = 1; % electron temperature in eV
q = 1.602e-19;
qm = q/9.11e-31;
epsi = 1/8.85e-12;
m = 9.11e-31;
n_0 = 1e20;
v = linspace(-1e7,1e7,Nv);
dv = v(2)-v(1);
dx = 1/(Nx-1);
dt = 5e-40;
M = sqrt(m/(2*3.1415*T*q))*exp(-m.*v.*v/(2*T*q));
mu1 = 0.3e7;
M_L = sqrt(m/(2*3.1415*T*q))*exp(-m.*(v-mu1).^2/(2*T*q));
mu2 = -0.3e7;
M_R = sqrt(m/(2*3.1415*T*q))*exp(-m.*(v-mu2).^2/(2*T*q));
E=1;
f = zeros(Nx,Nv);
f0 = zeros(Nx,Nv);
f1 = zeros(Nx,Nv);
fs = zeros(Nx,Nv,Nt);
for i=1:Nx
    f(i,:) = n_0/2*M_L*1.01+n_0/2*M_R;
    f0(i,:) =n_0*M;
    if mod(i,20) < 4 
      f(i,:) = n_0/2*M_L*0.96+n_0/2*M_R; 
    end
end
fs(:,:,1) = f;

for i=2:Nt
    f1 = zeros(Nx,Nv);
    
erho = zeros(Nx,1);
ionrho = zeros(Nx,1);

for j=1:Nx
   erho(j) = sum(f(j,:))*dv;
   ionrho(j) = sum(f0(j,:))*dv;
end

    erho = -erho*q/(dx);
    

        ionrho = ionrho*q/(dx);

        
        rho = erho+ionrho;
        
  % FIELDS - E field solver

A = zeros(Nx,Nx);

A(1,1) = -2;
A(1,2) = 1;



for l=2:Nx-1
    A(l,l-1) = 1;
    A(l,l) = -2;
    A(l,l+1) = 1;
end


A(Nx,Nx-1) = 1;
A(Nx,Nx) = -2;

phi = inv(A)*(-rho)*epsi*(dx*dx);

E = zeros(Nx,1);

E(1) = (phi(1)-phi(2))/(dx);

for j=2:Nx-1
    E(j) = (phi(j-1)-phi(j+1))/(2*dx);
end

E(Nx) = (phi(Nx-1)-phi(Nx))/(dx);
E;

A = zeros(Nx,Nv);
B = zeros(Nx,Nv);
A(1,:) = v.*(f(2,:) - f(Nx,:))/(2*dx);
B(:,1) = -qm*E.*(f(:,2)-f(:,Nv))/(2*dv);


for j=2:Nx-1
A(j,:) = v.*(f(j+1,:) - f(j-1,:))/(2*dx);
end
for k=2:Nv-1
B(:,k) = -qm*E.*(f(:,k+1)-f(:,k-1))/(2*dv);        
end

A(Nx,:) = v.*(f(1,:) - f(Nx-1,:))/(2*dx);
B(:,Nv) = -qm*E.*(f(:,1)-f(:,Nv))/(2*dv);

fsum = zeros(Nx,Nv);
fsum(1,1) = f(2,1)+f(Nx,1)+f(1,2)+f(1,Nv);
fsum(1,Nv) = f(2,Nv)+f(Nx,Nv)+f(1,1)+f(1,Nv-1);
fsum(1,2:Nv-1) = f(2,2:Nv-1)+f(Nx,2:Nv-1)+f(1,3:Nv)+f(1,1:Nv-2);
fsum(Nx,2:Nv-1) = f(1,2:Nv-1)+f(Nx-1,2:Nv-1)+f(Nx,3:Nv)+f(Nx,1:Nv-2);
fsum(2:Nx-1,1) = f(2:Nx-1,2)+f(2:Nx-1,Nv)+f(3:Nx,1)+f(1:Nx-2,1);
fsum(2:Nx-1,Nv) = f(2:Nx-1,2)+f(2:Nx-1,Nv-1)+f(3:Nx,Nv) + f(1:Nx-2,Nv);
fsum(Nx,1) = f(1,1)+f(Nx-1,1)+f(Nx,2)+f(Nx,Nv);
fsum(Nx,Nv) = f(1,Nv)+f(Nx-1,Nv)+f(Nx,1)+f(Nx,Nv-1);


for j=2:Nx-1
    for k=2:Nv-1
fsum(j,k) = f(j+1,k) +f(j-1,k)+f(j,k+1)+f(j,k-1);
    end
end
f1 = 1/4*fsum - (A+B)*dt;
fs(:,:,i) = f1;
f = f1;

 contourf(gridx,v,transpose(fs(:,:,i)),16)
 axis([0 1 -1e7 1e7])
 colormap(winter)
F(i-1) = getframe(gcf);
 i
end

%for n=1:Nt
%    plot(v,fs(1,:,n))
%    hold on
%end
numtimes=1;

fps=1;
figure
movie(F,numtimes,fps);

movie2avi(F, 'two_stream_LF.avi');
