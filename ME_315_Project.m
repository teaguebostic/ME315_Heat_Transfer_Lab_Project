%Teague Bostic
%ME 315 - Heat Transfer Lab
%Numerical Project
%Started: 11/6/2019
%Last Updated: 11/21/2019

clc; clear all; close all

%% define number of nodes

n = str2double(inputdlg('Number of Nodes: ','Nodes')); 

%% define parameters

k = 11.7; %W/(mK)

%% define dimensions of profile
x_mm = 75; %mm
y_mm = 25; %mm

x_m = x_mm*(10^-3); %m
y_m = y_mm*(10^-3); %m 

r = x_m/y_m; %nondimesional parameter

%% determine delta x and y
delx = x_m/(n-1);
dely = y_m/(n-1);

%% develop profile grid
x = linspace(0,x_m,n);
y = linspace(0,y_m,n);

[X,Y] = meshgrid(x,y);

%% develop temperature matrix
T = ones(n,n);
T = T*750;
%% define boundary conditions
T_LE = 800; %Celcius
h_LE = 200; %W/((m^2)K)
Bi_LE = (h_LE*dely)/k; %nondimensional

T_TE = 800; %Celcius
h_TE = 100; %W/((m^2)K)
Bi_TE = (h_TE*dely)/k; %nondimensional

T_top = 800; %Celcius
h_top = 150; %W/((m^2)K)
Bi_top = (h_top*delx)/k; %nondimensional

T_bot = 800; %Celcius
h_bot = 150; %W/((m^2)K)
Bi_bot = (h_bot*delx)/k; %nondimensional

%% determine energy sink

%Define Given Parameters
An = -15*(10^5);
theta = pi()/2; %rad

xn_mm = [12,28,43,55,67,67];
yn_mm = [0,0,0,0,5,-5]+12.5;

x_n = xn_mm*(10^-3);
y_n = yn_mm*(10^-3);

sigma_xn = [5,4.5,4,4,3,3]./1000;
sigma_yn = [3,2.5,2,2,1,1]./1000;

an = ((cos(theta)^2)./(2*((sigma_xn).^2)))+((sin(theta)^2)./(2*((sigma_yn).^2)));
bn = (-(sin(2*theta))./(4*((sigma_xn).^2)))+((sin(2*theta))./(4*((sigma_yn).^2)));
cn = ((sin(theta)^2)./(2*((sigma_xn).^2)))+((cos(theta)^2)./(2*((sigma_yn).^2)));

%Create an anyomous function of energy generation
e = @(X,Y,i) An.*exp(-((an(i).*((X-x_n(i)).^2))+(2.*bn(i).*((X-x_n(i))).*(Y-y_n(i)))+(cn(i).*((Y-y_n(i)).^2))));

%Solve for Enery Generation
e_gen = e(X,Y,1);
i = 2;
while i <= length(x_n)
e_gen = e_gen + e(X,Y,i);
i = i+1;
end

%Plot Energy Generation
figure(1)
contourf(X,Y,e_gen)
%% Define Enery Sink Matrix

e_sink = zeros(n,n);
%% Interior Node

for s = 2:1:n-1
    for d = 2:1:n-1
        xs = linspace(x(s)-(delx/2),x(s)+(delx/2),5);
        ys = linspace(y(d)-(dely/2),y(d)+(dely/2),5);
        [Xs, Ys] = meshgrid(xs,ys);
        
        es = e(Xs,Ys,1);
        
        i = 2;
        while i <= length(x_n)
            es = es + e(Xs,Ys,i);
            i = i+1;
        end
        I = trapz(ys,trapz(xs,es,2),1);
        V = delx*dely;
        e_sink(d,s) = I/V;
    end
end
%% Corner Node TL

s = 1;
d = 1;
xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
[Xs, Ys] = meshgrid(xs,ys);

es = e(Xs,Ys,1);

i = 2;
while i <= length(x_n)
    es = es + e(Xs,Ys,i);
    i = i+1;
end
I = trapz(ys,trapz(xs,es,2),1);
V = delx*dely;
e_sink(d,s) = I/V;
%% Corner Node TR

s = n;
d = 1;
xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
[Xs, Ys] = meshgrid(xs,ys);

es = e(Xs,Ys,1);

i = 2;
while i <= length(x_n)
    es = es + e(Xs,Ys,i);
    i = i+1;
end
I = trapz(ys,trapz(xs,es,2),1);
V = delx*dely;
e_sink(d,s) = I/V;
%% Corner Node BL

s = 1;
d = n;
xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
[Xs, Ys] = meshgrid(xs,ys);

es = e(Xs,Ys,1);

i = 2;
while i <= length(x_n)
    es = es + e(Xs,Ys,i);
    i = i+1;
end
I = trapz(ys,trapz(xs,es,2),1);
V = delx*dely;
e_sink(d,s) = I/V;
%% Corner Node BR

s = n;
d = n;
xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
[Xs, Ys] = meshgrid(xs,ys);

es = e(Xs,Ys,1);

i = 2;
while i <= length(x_n)
    es = es + e(Xs,Ys,i);
    i = i+1;
end
I = trapz(ys,trapz(xs,es,2),1);
V = delx*dely;
e_sink(d,s) = I/V;
%% Side Node LE

for d = 2:1:n-1
    s = 1;
    xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
    ys = linspace(y(d)-(dely/2),y(d)+(dely/2),5);
    [Xs, Ys] = meshgrid(xs,ys);
    
    es = e(Xs,Ys,1);
    
    i = 2;
    while i <= length(x_n)
        es = es + e(Xs,Ys,i);
        i = i+1;
    end
    I = trapz(ys,trapz(xs,es,2),1);
    V = delx*dely;
    e_sink(d,s) = I/V;
end
%% Side Node TE

for d = 2:1:n-1
    s = n;
    xs = linspace(x(s)-(delx/4),x(s)+(delx/4),5);
    ys = linspace(y(d)-(dely/2),y(d)+(dely/2),5);
    [Xs, Ys] = meshgrid(xs,ys);
    
    es = e(Xs,Ys,1);
    
    i = 2;
    while i <= length(x_n)
        es = es + e(Xs,Ys,i);
        i = i+1;
    end
    I = trapz(ys,trapz(xs,es,2),1);
    V = delx*dely;
    e_sink(d,s) = I/V;
end
%% Side Node Top

for s = 2:1:n-1
    d = 1;
    xs = linspace(x(s)-(delx/2),x(s)+(delx/2),5);
    ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
    [Xs, Ys] = meshgrid(xs,ys);
    
    es = e(Xs,Ys,1);
    
    i = 2;
    while i <= length(x_n)
        es = es + e(Xs,Ys,i);
        i = i+1;
    end
    I = trapz(ys,trapz(xs,es,2),1);
    V = delx*dely;
    e_sink(d,s) = I/V;
end
%% Side Node Bot

for s = 2:1:n-1
    d = n;
    xs = linspace(x(s)-(delx/2),x(s)+(delx/2),5);
    ys = linspace(y(d)-(dely/4),y(d)+(dely/4),5);
    [Xs, Ys] = meshgrid(xs,ys);
    
    es = e(Xs,Ys,1);
    
    i = 2;
    while i <= length(x_n)
        es = es + e(Xs,Ys,i);
        i = i+1;
    end
    I = trapz(ys,trapz(xs,es));
    V = delx*dely;
    e_sink(d,s) = I/V;
end

figure(2)
contourf(X,Y,e_sink)
%% solve nodal equations
z = 1;

while 1
    T_old = T;
    
    for i = 2:1:(n-1)
        for j = 2:1:(n-1)
            T(i,j) = (T(i,j-1)+((r^2)*T(i+1,j))+T(i,j+1)+((r^2)*T(i-1,j))+((r)*((e_sink(i,j)*delx*dely)/k)))/(2*(1+(r^2)));
        end
    end
    
    for j = 1
        i = 1;
        T(i,j) = ((r*Bi_LE*T_LE)+((r^2)*T(i+1,j))+(T(i,j+1))+((r*Bi_top*T_top)+((r)*((e_sink(i,j)*delx*dely)/(2*k)))))/((r*Bi_LE)+(r^2)+1+(r*Bi_top));
    end
    
    for j = n
        i = 1;
        T(1,j) = ((r*Bi_TE*T_TE)+((r^2)*T(i+1,j))+T(i,j-1)+(r*Bi_top*T_top)+((r)*((e_sink(i,j)*delx*dely)/(2*k))))/((r*Bi_TE)+(r^2)+1+(r*Bi_top));
    end
    
    for j = 1
        i = n;
        T(i,j) = ((r*Bi_LE*T_LE)+((r^2)*T(i-1,j))+T(i,j+1)+(r*Bi_bot*T_bot)+((r)*((e_sink(i,j)*delx*dely)/(2*k))))/((r*Bi_LE)+(r^2)+1+(r*Bi_bot));
    end
    
    for j = n
        i = n;
        T(i,j) = ((r*Bi_TE*T_TE)+((r^2)*T(i-1,j))+T(i,j-1)+(r*Bi_bot*T_bot)+((r)*((e_sink(i,j)*delx*dely)/(2*k))))/((r*Bi_TE)+(r^2)+1+(r*Bi_bot));
    end
    
    for i = 2:1:(n-1)
        j = 1;
        T(i,j) = ((2*r*Bi_LE*T_LE)+(r^2)*T(i+1,j)+(2*T(i,j+1))+(r^2)*T(i-1,j)+((r)*((e_sink(i,j)*delx*dely)/k)))/((2*r*Bi_LE)+(2*(r^2))+2);
    end
    
    for i = 2:1:(n-1)
        j = n;
        T(i,j) = ((2*r*Bi_TE*T_TE)+(r^2)*T(i+1,j)+(2*T(i,j-1))+(r^2)*T(i-1,j)+((r)*((e_sink(i,j)*delx*dely)/k)))/((2*r*Bi_TE)+(2*(r^2))+2);
    end
    
    for j = 2:1:(n-1)
        i = 1;
        T(i,j) = ((2*r*Bi_top*T_top)+(r^2)*T(i,j+1)+(2*T(i+1,j))+(r^2)*T(i,j-1)+((r)*((e_sink(i,j)*delx*dely)/k)))/((2*r*Bi_top)+(2*(r^2))+2);
    end
    
    for j = 2:1:(n-1)
        i = n;
        T(i,j) = ((2*r*Bi_bot*T_bot)+(r^2)*T(i,j+1)+(2*T(i-1,j))+(r^2)*T(i,j-1)+((r)*((e_sink(i,j)*delx*dely)/k)))/((2*r*Bi_bot)+(2*(r^2))+2);
    end
    
    T_diff = (T(n,n) - T_old(n,n))/T_old(n,n);
    z = z+1;
    
    if T_diff < 0.000000001
        break
    end
end


%% Contour and Mesh plots for Temperature

figure(3)
mesh(X,Y,T);
zlim([725 825])

figure(4)
contourf(X,Y,T)