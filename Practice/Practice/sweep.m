
%given data

sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;

%spatial discretization

X=4;
Y=4;

dx=0.05;
dy=0.05;

x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);

%angular discretization
N=16;

[mu, w]=lgwt(N,-1,1);

cell_phi=cell(N,1);
eta=cell(N,1);
del_phi=zeros(N,1);
division_per_mu_direction=0;
for i=1:N/2
    division_per_mu_direction=division_per_mu_direction+4;
    del_phi(i,1)=2*pi/division_per_mu_direction;
    cell_phi{i,1}=(del_phi(i,1)/2:del_phi(i,1):2*pi);
    eta{i,1}=sqrt(1-mu(i,1)*mu(i,1))*cos(cell_phi{i,1});
    
end

for i=N/2+1:N
    del_phi(i,1)=2*pi/division_per_mu_direction;
    cell_phi{i,1}=(del_phi(i,1)/2:del_phi(i,1):2*pi);
    division_per_mu_direction=division_per_mu_direction-4;
    eta{i,1}=sqrt(1-mu(i,1)*mu(i,1))*cos(cell_phi{i,1});
end

psi_out=zeros(n_x,n_y);
scaler_flux=zeros(n_x-1,n_y-1);

for i=1:n_y
    for j=1:n_x
        for p=1:N/2
            for q=length(eta{p})



    