% FIRST NAME, LAST NAME = SELİN, TÜRKMEN

clear all
close all

%% Question 1
e0 = 8.854*(10^-12);
c = 1/(4*pi*e0);
q1 = 3*10^-10;
q2 = 6*10^-10;

[x,y] = meshgrid(-0.75:0.01:0.75,-0.75:0.01:0.75);
V = c.*((q1./sqrt((x-0.5).^2+y.^2))+(q2./sqrt((x+0.5).^2+y.^2)));

figure
surf(x,y,V),xlabel('x'),ylabel('y'),zlabel('V'), title('Electric Potential Field')
figure
meshc(x,y,V),xlabel('x'),ylabel('y'),zlabel('V'), title('Electric Potential Field')

%% Question 2
h = animatedline;
x = linspace(-100,350);
xlabel('Distance(m)')
ylabel('Amplitude(V)')
title('Free-Space EM Wave Propagation')
legend('Amplitude= -2.4317 V')

for k = 1:length(x)
    y = cos(2*pi*5*x(k));
    addpoints(h,x(k),y);
    drawnow limitrate
end
drawnow
%% Question 3

han = animatedline;

m = 1.5;
h0 = 26;
v0 = 18;
g = 9.8;
E(1) = m*g*h0 + m*v0^2/2;
hmax(1) = E(1)/(m*g);
t2(1) = v0/g + sqrt(2*hmax(1)/g);
t0= t2(1);
t = [0:0.01:t0];
h = h0 + v0*t -(g*t.^2)/2;

addpoints(han,t,h);
drawnow limitrate

i=1;
while E(i)>10
    E(i+1)= 0.8*E(i);
    hmax(i+1) = E(i+1)/(m*g);
    v(i)= sqrt(2*E(i+1)/m);
    t2(i+1) = t2(i)+2*sqrt(2*hmax(i+1)/g);
    tr2=t2(i+1);
    tr1=t2(i);
    t= [tr1:0.01:tr2];
    h1 = (v(i)*(t-t2(i)))-((g*(t-t2(i)).^2)/2);
    
hold on
addpoints(han,t,h1);
drawnow limitrate

i=i+1;
end
drawnow


%% HW END