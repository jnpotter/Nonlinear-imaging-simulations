% nonlinear imaging field simulation
clear

global mu rho K dx dy n m A B C gpu_on
global duxdxh duxdyh duydxh duydyh dduxdxdxh dduxdxdyh dduxdydyh dduydydyh dduydxdyh dduydxdxh
global duxdx duxdy duydx duydy dduxdxdx dduxdxdy dduxdydy dduydydy dduydxdy dduydxdx

gpu_on=1;

%material properties

E=70e9;
v=.35;

K=E/3/(1-2*v); %bulk modulus
mu=E/2/(1+v); %shear modulus
rho=2700;

A=-100e9;B=-100e9;C=-100e9; %TOECS

cl=sqrt((K+4/3*mu)/rho);
ct=sqrt(mu/rho);

array_fname = 'Imasonic 1D 64els 5.00MHz 0.63mm pitch.mat';

%test parameters
T=5*1e-6; %test length
t0=0e-6;
dt=1/5e6/25; %time step

dx=cl/5e6/25; %spatial step
dz=dx;
dy=dx;

xrange=[-10 10]*1e-3;
zrange=[-2 30]*1e-3; %15

%focal point
xf=0; zf= 1.5e-2;

%input signal
centre_freq = 5e6;
wh=2/3*centre_freq; %half power frequency
delta=10e-9; %excitation amplitude

%timebase details (for pulse)
time_step = 1/1000e6;
time_pts =100000;% round(T/time_step);

alpha=-log(0.5)/(centre_freq-wh)^2; %gaussian coefficient
freq=((1:time_pts)/time_pts/time_step)';
time=time_step:time_step:time_pts*time_step;

in_freq_spec=exp(-alpha*(freq-centre_freq).^2).*exp(-i*2*pi*freq*.5e-6); %gaussian input pulse
in_time_sig=imag(ifft(in_freq_spec));
in_freq_spec=in_freq_spec/max(abs(hilbert(in_time_sig)))*delta;
in_time_sig=imag(ifft(in_freq_spec));

x_arr=xrange(1):dx:xrange(2);
z_arr=zrange(1):dz:zrange(2);
[x z]=meshgrid(x_arr,z_arr);
xn=length(x_arr);zn=length(z_arr);

zero_els=find(z<0);

tmp = load(array_fname);
exp_data.array = tmp.array;

exp_data.num_els = length(exp_data.array.el_xc);
a= exp_data.array.el_x1(1)- exp_data.array.el_x2(1); %element width

%focal delay arr
delay_arr=sqrt((exp_data.array.el_xc-xf).^2 + (exp_data.array.el_zc-zf).^2)/cl;
delay_arr=delay_arr-min(delay_arr);

%max directivity value
theta_mat=0;
D_mat=sin(pi*a*sin(theta_mat)/(cl/centre_freq))./(pi*a*sin(theta_mat)/(cl/centre_freq));
zeta=sin(theta_mat);
F0_mat=(2*zeta.^2 -(cl/ct)^2).^2 - 4*zeta.^2.*(zeta.^2 -1).^(0.5) .*(zeta.^2 -(cl/ct)^2).^(0.5);
DL_max=((cl/ct)^2 -2*sin(theta_mat).^2).*cos(theta_mat)./F0_mat;

n=size(x,1);m=size(x,2);

if gpu_on==1 %memory preallocation
    
    duxdxh=gpuArray(zeros(n,m)); duxdyh=gpuArray(zeros(n,m));  duydxh=gpuArray(zeros(n,m));   duydyh=gpuArray(zeros(n,m));  dduxdxdxh=gpuArray(zeros(n,m));   dduxdxdyh=gpuArray(zeros(n,m));   dduxdydyh=gpuArray(zeros(n,m));    dduydydyh=gpuArray(zeros(n,m));    dduydxdyh=gpuArray(zeros(n,m));    dduydxdxh=gpuArray(zeros(n,m));
    duxdx=gpuArray(zeros(n,m)); duydy=gpuArray(zeros(n,m)); duxdy=gpuArray(zeros(n,m)); duydx=gpuArray(zeros(n,m)); dduxdxdx=gpuArray(zeros(n,m)); dduxdxdy=gpuArray(zeros(n,m));dduxdydy=gpuArray(zeros(n,m));dduydydy=gpuArray(zeros(n,m));dduydxdy=gpuArray(zeros(n,m));   dduydxdx=gpuArray(zeros(n,m));
    
    ux=gpuArray(zeros(size(x)));
    r_mat=gpuArray(zeros(size(x)));
    theta_mat=gpuArray(zeros(size(x)));

    in_freq_spec=gpuArray(in_freq_spec);
    freq=gpuArray(freq);
    
    sin_mat=gpuArray(zeros(size(x)));
    cos_mat=gpuArray(zeros(size(x)));
    
    in_time_sig_prop=gpuArray(zeros(size(in_time_sig)));
    delay_arr=gpuArray(delay_arr);
    
    uy=gpuArray(ux);ux_dot=gpuArray(ux);uy_dot=gpuArray(ux);Fx=gpuArray(ux);Fy=gpuArray(ux); psf=ux;%Initial states
    
else
    
    duxdxh=(zeros(n,m));   duxdyh=(zeros(n,m));    duydxh=(zeros(n,m));  duydyh=(zeros(n,m));   dduxdxdxh=(zeros(n,m));   dduxdxdyh=(zeros(n,m));   dduxdydyh=(zeros(n,m));   dduydydyh=(zeros(n,m)); dduydxdyh=(zeros(n,m));   dduydxdxh=(zeros(n,m));
    duxdx=zeros(n,m); duydy=zeros(n,m); duxdy=zeros(n,m); duydx=zeros(n,m);  dduxdxdx=zeros(n,m);dduxdxdy=zeros(n,m); dduxdydy=zeros(n,m);  dduydydy=zeros(n,m);  dduydxdy=zeros(n,m);    dduydxdx=zeros(n,m);
    ux=zeros(size(x));uy=ux;ux_dot=ux;uy_dot=ux;Fx=ux;Fy=ux; psf=ux;%Initial states
    
end

figure(1);subplot(1,2,1);fig_han1=imagesc(x_arr,z_arr,(ux));axis equal;axis tight;colorbar;title('Linear field');%caxis([0 3e-8]);
figure(1);subplot(1,2,2);fig_han2=imagesc(x_arr,z_arr,(uy));colorbar;axis equal;axis tight;colorbar;title('Nonlinear field');%caxis([0 3.5e-11])

b_length=0.005; b_weight=10000;
left_els=find(x<xrange(1)+b_length);
right_els=find(x>xrange(2)-b_length);
top_els=find(z<zrange(1)+b_length);
bottom_els=find(z>zrange(2)-b_length);

left_mat=exp(-b_weight*(x(left_els)-(xrange(1)+b_length)).^2);
right_mat=exp(-b_weight*(x(right_els)-(xrange(2)-b_length)).^2);
top_mat=exp(-b_weight*(z(top_els)-(zrange(1)+b_length)).^2);
bottom_mat=exp(-b_weight*(z(bottom_els)-(zrange(2)-b_length)).^2);

ii=0;

for t=t0:dt:T; %time point
    
    tic;
    ii=ii+1;
    
    ux_arr=zeros(size(x));
    uz_arr=zeros(size(x));
    
    for N=1:exp_data.num_els %evaluating linear field
         
        theta_mat=atan((x-exp_data.array.el_xc(N))./(z-exp_data.array.el_zc(N)));
        r_mat=sqrt((x-exp_data.array.el_xc(N)).^2+(z-exp_data.array.el_zc(N)).^2);
        
        D_mat=sin(pi*a*sin(theta_mat)/(cl/centre_freq))./(pi*a*sin(theta_mat)/(cl/centre_freq)); %directivity
        
        zeta=sin(theta_mat);
        F0_mat=(2*zeta.^2 -(cl/ct)^2).^2 - 4*zeta.^2.*(zeta.^2 -1).^(0.5) .*(zeta.^2 -(cl/ct)^2).^(0.5);
        
        DL_mat=((cl/ct)^2 -2*sin(theta_mat).^2).*cos(theta_mat)./F0_mat;
        DL_mat=DL_mat/DL_max;
        
        in_time_sig_prop=imag(ifft(in_freq_spec.*exp(-1i*2*pi*freq*(t+delay_arr(N)))));   %propagating pulse
        u_arr=interp1(time, in_time_sig_prop, r_mat/cl)./sqrt(r_mat)*(sqrt(cl/centre_freq)).*D_mat.*DL_mat; %divergence and directivity
        u_arr(zero_els)=-u_arr(zero_els);
        
        ux_arr=ux_arr+u_arr.*sin(theta_mat);
        uz_arr=uz_arr+u_arr.*cos(theta_mat);
        
    end
    
    %crude absorbing boundaries
    ux(left_els)=ux(left_els).*left_mat;
    ux(right_els)=ux(right_els).*right_mat;
    ux(top_els)=ux(top_els).*top_mat;
    ux(bottom_els)=ux(bottom_els).*bottom_mat;
    
    uy(left_els)=uy(left_els).*left_mat;
    uy(right_els)=uy(right_els).*right_mat;
    uy(top_els)=uy(top_els).*top_mat;
    uy(bottom_els)=uy(bottom_els).*bottom_mat;
    
    %%Nonlinear time-stepping
    [Fx, Fy]=non_linear_forcing_computation1(ux_arr, uz_arr); %computes nonlinear forcing at current time-step
    
    %%4th order 2D runge-kutta time-stepping routine
    [ux_dotdot, uy_dotdot]=nonlin_wave_eqs_high_order(ux, uy, Fx, Fy, ux_dot, uy_dot);
    k11_ux=ux_dot*dt;
    k12_ux=ux_dotdot*dt;
    k11_uy=uy_dot*dt;
    k12_uy=uy_dotdot*dt;
    
    [ux_dotdot, uy_dotdot]=nonlin_wave_eqs_high_order(ux+ k11_ux/2, uy+k11_uy/2, Fx, Fy, ux_dot+k12_ux/2, uy_dot+k12_uy/2);
    k21_ux=(ux_dot +k12_ux/2)*dt;
    k22_ux=ux_dotdot*dt;
    k21_uy=(uy_dot+k12_uy/2)*dt;
    k22_uy=uy_dotdot*dt;
    
    [ux_dotdot, uy_dotdot]=nonlin_wave_eqs_high_order(ux+ k21_ux/2, uy+k21_uy/2, Fx, Fy, ux_dot+k22_ux/2, uy_dot+k22_uy/2);
    k31_ux=(ux_dot +k22_ux/2)*dt;
    k32_ux=ux_dotdot*dt;
    k31_uy=(uy_dot+k22_uy/2)*dt;
    k32_uy=uy_dotdot*dt;
    
    [ux_dotdot, uy_dotdot]=nonlin_wave_eqs_high_order(ux+ k31_ux, uy+k31_uy, Fx, Fy, ux_dot+k32_ux, uy_dot+k32_uy);
    k41_ux=(ux_dot +k32_ux)*dt;
    k42_ux=ux_dotdot*dt;
    k41_uy=(uy_dot+k32_uy)*dt;
    k42_uy=uy_dotdot*dt;
    
    ux_dot=ux_dot+(k12_ux+2*k22_ux +2*k32_ux + k42_ux)/6;
    uy_dot=uy_dot+(k12_uy+2*k22_uy +2*k32_uy + k42_uy)/6;
    
    ux=ux+(k11_ux+2*k21_ux +2*k31_ux + k41_ux)/6; %displacements of nonlinear field
    uy=uy+(k11_uy+2*k21_uy +2*k31_uy + k41_uy)/6;
    
    ux(zero_els)=0;uy(zero_els)=0; %remove field behind array
    ux_arr(zero_els)=0;uz_arr(zero_els)=0; %remove field behind array
    
    
    psf=psf+(Fx.*ux_dot + Fy.*uy_dot)*dt;
    
    set(fig_han1,'CData', gather(sqrt(uz_arr.^2+ux_arr.^2)));
    set(fig_han2,'CData', gather(sqrt(uy.^2+ux.^2)));
    drawnow
    
    %  nonlinear_mov(ii)=getframe;
    
    disp(toc)
    
end

figure;imagesc(x_arr,z_arr, psf)




