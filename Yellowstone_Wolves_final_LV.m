%%Impact of Wolf Reintroduction on Yellowstone Elk Population Model
%%1995-Present%%
clc
clear all
close all

%%Initial Conditions%%
niter = 10000;  
dt = 0.001;    

A(1) = 1.6899e4; % initial elk population (prey)
B(1) = 31; % initial wolf population (predator)
r1 = 1.4; %elk population growth rate (accounts for non-predator death)
r2 = 0.04855; %predation rate
r3 = 0.00009511; %predator rate of birth per per prey
r4 = 0.2411; %predator rate of death

%%2-step Runge-Kutta%%
for i=1:niter
    time(i+1)=i*dt;
    
    %Population cannot be negative
    A(A < 0) = 0;
    B(B < 0) = 0;
   
    K1A=fa(A(i),B(i),r1,r2,dt);
    K1B=fb(A(i),B(i),r3,r4,dt);
   
    K2A=fa(A(i)+0.5*K1A,B(i)+0.5*K1B,r1,r2,dt);
    K2B=fb(A(i)+0.5*K1A,B(i)+0.5*K1B,r3,r4,dt);
    
    K3A=fa(A(i)+0.5*K2A,B(i)+0.5*K2B,r1,r2,dt);
    K3B=fb(A(i)+0.5*K2A,B(i)+0.5*K2B,r3,r4,dt);
    
    K4A=fa(A(i)+0.5*K3A,B(i)+0.5*K3B,r1,r2,dt);
    K4B=fb(A(i)+0.5*K3A,B(i)+0.5*K3B,r3,r4,dt);
   
    A(i+1)=A(i)+(K1A+(2*K2A)+(2*K3A)+K4A)/6; 
    B(i+1)=B(i)+(K1B+(2*K2B)+(2*K3B)+K4B)/6;
       
end

disp (max(A));
disp (max(B));

%%plots%%
subplot(2,1,1);
plot(time,A)
title(['Yellowstone Elk Population Projection with Hunting from 2015'], 'fontsize', 10)
xlabel(['Time (years)'], 'fontsize', 8)
ylabel(['Elk Population.'], 'fontsize', 8)


subplot(2,1,2);
plot(time,B)
title(['Yellowstone Wolf Population Projection with Hunting from 2015'], 'fontsize', 10)
xlabel(['Time (years)'], 'fontsize', 8)
ylabel(['Wolf Population'], 'fontsize', 8)

%%Runge Kutta functions%%
%Pop. A
function KA=fa(A,B,r1,r2,dt)
KA=dt*((r1*A)-(r2*A*B));
end

%Pop. B
function KB=fb(A,B,r3,r4,dt)
KB=dt*((r3*A*B)-(r4*B));
end 


