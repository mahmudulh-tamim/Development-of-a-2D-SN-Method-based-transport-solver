
function [k_new, flux_new]=power_iteration()
%given data
tol=10^(-7);


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
tot_angular_direction_count=N*(N+2)/2;

%%

flux_old=(1/sum(sum(ones(n_x-1,n_y-1)*dx*dy*nu_sigma_f)))*ones(n_x-1,n_y-1);
k_old=1;

flux_new_half=source_iteration(flux_old,k_old);
k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

flux_new=k_old/k_new*flux_new_half;

iteration=1;

while abs(k_new-k_old)>tol
    flux_old=flux_new;
    k_old=k_new;
    flux_new_half=source_iteration(flux_old,k_old);
    k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

    flux_new=k_old/k_new*flux_new_half;

    iteration=1+iteration;
end
mesh_centre_x=(0.05/2:0.05:4)';
mesh_centre_y=(0.05/2:0.05:4)';

mesh(mesh_centre_x, mesh_centre_y,flux_new);

iteration