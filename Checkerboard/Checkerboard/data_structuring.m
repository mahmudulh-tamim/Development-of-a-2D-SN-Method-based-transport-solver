%given data

%first type of cell
sigma_t_1=1;
sigma_s_1=0.7;
nu_sigma_f_1=0.24;

%second type of cell
sigma_t_2=1;
sigma_s_2=0.7;
nu_sigma_f_2=0.39;

%spatial discretization

X=(0:4:24);
Y=(0:4:24);

dx=0.2;
dy=0.2;

x_representative_of_single_cell=(0:dx:4)';
y_representative_of_single_cell=(0:dy:4)';
number_of_xmesh_of_a_single_cell=length(x_representative_of_single_cell)-1;
number_of_ymesh_of_a_single_cell=length(y_representative_of_single_cell)-1;
x=(0:dx:24)';
y=(0:dy:24)';


number_of_xmesh=length(x)-1;
number_of_ymesh=length(y)-1;

%initialization of the material geometry
sigma_t=zeros(number_of_xmesh,number_of_ymesh);
sigma_s=zeros(number_of_xmesh,number_of_ymesh);
nu_sigma_f=zeros(number_of_xmesh, number_of_ymesh);

flag1=true;
flag2=true;
start_index_1=1;
end_index_1=number_of_ymesh_of_a_single_cell;
start_index_2=1;
end_index_2=number_of_xmesh_of_a_single_cell; % doesn't matter as number_of_ymesh_of_a_single_cell=number_of_xmesh_of_a_single_cell
for i=1:6
    if flag1==true
        for j=1:6
            if flag2==true
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_1;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_1;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_1;
                flag2=false;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            elseif flag2==false
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_2;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_2;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_2;
                flag2=true;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            end
        end
        flag1=false;
        start_index_1=end_index_1+1;
        end_index_1=end_index_1+number_of_ymesh_of_a_single_cell;
        start_index_2=1;
        end_index_2=number_of_xmesh_of_a_single_cell;

    elseif flag1==false
        for j=1:6
            if flag2==true
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_2;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_2;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_2;
                flag2=false;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            elseif flag2==false
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_1;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_1;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_1;
                flag2=true;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            end
        end
        flag1=true;
        start_index_1=end_index_1+1;
        end_index_1=end_index_1+number_of_ymesh_of_a_single_cell;
        start_index_2=1;
        end_index_2=number_of_xmesh_of_a_single_cell;
    end
end




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
