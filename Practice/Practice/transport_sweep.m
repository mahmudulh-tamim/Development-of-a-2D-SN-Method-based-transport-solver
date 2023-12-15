function scaler_flux=transport_sweep(Q)
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
tot_angular_direction_count=N*(N+2)/2;
mu=[0.138956875067780344591732;
    0.392289261444811712294197;
    0.537096561300879079878296;
    0.650426450628771770509703;
    0.746750573614681064580018;
    0.831996556910044145168291;
    0.909285500943725291652116;
    0.980500879011739882135849];
eta=[0.138956875067780344591732;
    0.392289261444811712294197;
    0.537096561300879079878296;
    0.650426450628771770509703;
    0.746750573614681064580018;
    0.831996556910044145168291;
    0.909285500943725291652116;
    0.980500879011739882135849];
w=[0.0489872391580385335008367;
   0.0413295978698440232405505;
   0.0203032007393652080748070;
   0.0265500757813498446015484;
   0.0379074407956004002099321;
   0.0135295047786756344371600;
   0.0326369372026850701318409;
   0.0103769578385399087825920];
w_index = [1
         2
         3
         4
         4
         3
         2
         1
         2
         5
         6
         7
         6
         5
         2
         3
         6
         8
         8
         6
         3
         4
         7
         8
         7
         4
         4
         6
         6
         4
         3
         5
         3
         2
         2
         1];

psi_cell=zeros(n_x-1,n_y-1,tot_angular_direction_count);
row_interface_psi=zeros(n_y-1,n_x,tot_angular_direction_count);
col_interface_psi=zeros(n_y,n_x-1,tot_angular_direction_count);
scaler_flux=zeros(n_x-1,n_y-1);
ang_direct_count=0;

%% sweep from left to right and bottom to top
for j=1:n_y
    for i=1:n_x
        ang_direct_count=1;
        for q=1:N/2
            for p=1:N/2
                if (p+q>=2 && p+q<10)
                    if(j==1)
                        
                        if (i==1)
                            row_interface_psi(j,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i-1,ang_direct_count)=0;
                            psi_cell(i-1,j,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j,i-1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i-1,ang_direct_count)+Q(i-1,j))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j,i,ang_direct_count)=2*psi_cell(i-1,j,ang_direct_count)-row_interface_psi(j,i-1,ang_direct_count);
                            scaler_flux(i-1,j)=0.25*w(w_index(ang_direct_count,1))*psi_cell(i-1,j,ang_direct_count)+scaler_flux(i-1,j);

                        end
                        
                    elseif (j~=n_y)
                        
                        if (i==1)
                            row_interface_psi(j,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i-1,ang_direct_count)=2*psi_cell(i-1,j-1,ang_direct_count)-col_interface_psi(j-1,i-1,ang_direct_count);
                            psi_cell(i-1,j,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j,i-1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i-1,ang_direct_count)+Q(i-1,j))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j,i,ang_direct_count)=2*psi_cell(i-1,j,ang_direct_count)-row_interface_psi(j,i-1,ang_direct_count);
                            scaler_flux(i-1,j)=0.25*w(w_index(ang_direct_count,1))*psi_cell(i-1,j,ang_direct_count)+scaler_flux(i-1,j);
                        end

                    else 
                        if (i~=1)
                            col_interface_psi(j,i-1,ang_direct_count)=2*psi_cell(i-1,j-1,ang_direct_count)-col_interface_psi(j-1,i-1,ang_direct_count);
                        end    
                    end
                        

                   ang_direct_count=ang_direct_count+1; 
                end
            end
        end
    end
end

%% sweep from right to left and bottom to top
for j=1:n_y
    for i=n_x:-1:1
        ang_direct_count=37;
        for q=1:N/2
            for p=1:N/2
                if (p+q>=2 && p+q<10)
                    if(j==1)
                        
                        if (i==n_x)
                            row_interface_psi(j,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i,ang_direct_count)=0;
                            psi_cell(i,j,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j,i+1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i,ang_direct_count)+Q(i,j))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j,ang_direct_count)-row_interface_psi(j,i+1,ang_direct_count);
                            scaler_flux(i,j)=0.25*w(w_index(ang_direct_count-36,1))*psi_cell(i,j,ang_direct_count)+scaler_flux(i,j);

                        end
                        
                    elseif (j~=n_y)
                        
                        if (i==n_x)
                            row_interface_psi(j,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j-1,ang_direct_count)-col_interface_psi(j-1,i,ang_direct_count);
                            psi_cell(i,j,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j,i+1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i,ang_direct_count)+Q(i,j))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j,ang_direct_count)-row_interface_psi(j,i+1,ang_direct_count);
                            scaler_flux(i,j)=0.25*w(w_index(ang_direct_count-36,1))*psi_cell(i,j,ang_direct_count)+scaler_flux(i,j);
                        end

                    else 
                        if (i~=n_x)
                            col_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j-1,ang_direct_count)-col_interface_psi(j-1,i,ang_direct_count);
                        end    
                    end
                        

                   ang_direct_count=ang_direct_count+1; 
                end
            end
        end
    end
end

%% sweep from left to right and top to bottom
for j=n_y:-1:1
    for i=1:n_x
        ang_direct_count=73;
        for q=1:N/2
            for p=1:N/2
                if (p+q>=2 && p+q<10)
                    if(j==n_y)
                        
                        if (i==1)
                            row_interface_psi(j-1,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i-1,ang_direct_count)=0;
                            psi_cell(i-1,j-1,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j-1,i-1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i-1,ang_direct_count)+Q(i-1,j-1))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j-1,i,ang_direct_count)=2*psi_cell(i-1,j-1,ang_direct_count)-row_interface_psi(j-1,i-1,ang_direct_count);
                            scaler_flux(i-1,j-1)=0.25*w(w_index(ang_direct_count-72,1))*psi_cell(i-1,j-1,ang_direct_count)+scaler_flux(i-1,j-1);

                        end
                        
                    elseif (j~=1)
                        
                        if (i==1)
                            row_interface_psi(j-1,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i-1,ang_direct_count)=2*psi_cell(i-1,j,ang_direct_count)-col_interface_psi(j+1,i-1,ang_direct_count);
                            psi_cell(i-1,j-1,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j-1,i-1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i-1,ang_direct_count)+Q(i-1,j-1))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j-1,i,ang_direct_count)=2*psi_cell(i-1,j-1,ang_direct_count)-row_interface_psi(j-1,i-1,ang_direct_count);
                            scaler_flux(i-1,j-1)=0.25*w(w_index(ang_direct_count-72,1))*psi_cell(i-1,j-1,ang_direct_count)+scaler_flux(i-1,j-1);
                        end

                    else 
                        if (i~=1)
                            col_interface_psi(j,i-1,ang_direct_count)=2*psi_cell(i-1,j,ang_direct_count)-col_interface_psi(j+1,i-1,ang_direct_count);
                        end    
                    end
                        

                   ang_direct_count=ang_direct_count+1; 
                end
            end
        end
    end
end

%% sweep from right to left and top to bottom
for j=n_y:-1:1
    for i=n_x:-1:1
        ang_direct_count=109;
        for q=1:N/2
            for p=1:N/2
                if (p+q>=2 && p+q<10)
                    if(j==n_y)
                        
                        if (i==n_x)
                            row_interface_psi(j-1,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i,ang_direct_count)=0;
                            psi_cell(i,j-1,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j-1,i+1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i,ang_direct_count)+Q(i,j-1))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j-1,i,ang_direct_count)=2*psi_cell(i,j-1,ang_direct_count)-row_interface_psi(j-1,i+1,ang_direct_count);
                            scaler_flux(i,j-1)=0.25*w(w_index(ang_direct_count-108,1))*psi_cell(i,j-1,ang_direct_count)+scaler_flux(i,j-1);

                        end
                        
                    elseif (j~=1)
                        
                        if (i==n_x)
                            row_interface_psi(j-1,i,ang_direct_count)=0;
                        else
                            col_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j,ang_direct_count)-col_interface_psi(j+1,i,ang_direct_count);
                            psi_cell(i,j-1,ang_direct_count)=(2*mu(p,1)/dx*row_interface_psi(j-1,i+1,ang_direct_count)+2*eta(q,1)/dy*col_interface_psi(j,i,ang_direct_count)+Q(i,j-1))/(sigma_t+2*mu(p,1)/dx+2*eta(q,1)/dy);
                            row_interface_psi(j-1,i,ang_direct_count)=2*psi_cell(i,j-1,ang_direct_count)-row_interface_psi(j-1,i+1,ang_direct_count);
                            scaler_flux(i,j-1)=0.25*w(w_index(ang_direct_count-108,1))*psi_cell(i,j-1,ang_direct_count)+scaler_flux(i,j-1);
                        end

                    else 
                        if (i~=n_x)
                            col_interface_psi(j,i,ang_direct_count)=2*psi_cell(i,j,ang_direct_count)-col_interface_psi(j+1,i,ang_direct_count);
                        end    
                    end
                        

                   ang_direct_count=ang_direct_count+1; 
                end
            end
        end
    end
end



