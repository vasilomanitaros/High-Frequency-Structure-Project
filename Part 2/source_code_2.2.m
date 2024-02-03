Ien=[1 -1 1 -1 1 -1 1 -1];
I=[1 1 1 1 1 1 1 1];
Ipick=Ien;
%Vertical
phi=0;
theta=[0:180];
ef=efield(deg2rad([0:180]),deg2rad([-180:180]),Ipick);
maxe=max(ef, [], 'all');
for i=1:length(theta)
    evert(i)=efield(deg2rad(theta(i)),phi,Ipick);
end
evertmax=max(evert);
figure(1)
polarplot(deg2rad(theta),evert/evertmax);
title('Vertical Polar Plot','0<theta<180, Z-XY Pane')

%Horizontal
phi=[0:360];
theta=pi/2;
for i=1:length(phi)
    ehor(i)=efield(theta,deg2rad(phi(i)),Ipick);
end
ehormax=max(ehor)
figure(2)
polarplot(deg2rad(phi),ehor/ehormax);
title('Horizontal Polar Plot','0<phi<360, X-Y Pane')

%3Dimensional
figure(4)
syms th ph
rho=efield(th,ph,Ipick);
x = rho*sin(th)*cos(ph);
y = rho*sin(th)*sin(ph);
z = rho*cos(th);
fsurf(x,y,z,[0 2*pi 0 pi])
xlabel('X')
ylabel('Y')
zlabel('Z')

figure(3)
patternCustom(ef.',[0:180],[0:360]);


%Efield Calculation
function efield=efield(theta,phi,I)
%Constants
f=10^9;
lambda=physconst('LightSpeed')/f;
k=2*pi/lambda;
d=3*lambda/4;
r=10;%does not matter
%Antenna Initialization
antennas=[-7*d/2 -5*d/2 -3*d/2 -d/2 +d/2 +3*d/2 +5*d/2 +7*d/2];
efieldfinal=zeros(length(theta),length(phi));
for i=1:8
    for ph=1:length(phi)
        for th=1:length(theta)
            if (theta(th)==0 || theta(th)==pi) efield(th,ph)=0; 
                else
            efield(th,ph)=j*60*I(i)/(r*sin(theta(th)))*exp(-j*k*(r-antennas(i)*cos(phi(ph)))*sin(theta(th)))*cos(pi/2*cos(theta(th)));
            end
        end
    end
    efieldfinal=efieldfinal+efield;
end
efield=abs(efieldfinal);
end


%for one antenna
% efield=[j.*(120*pi/(4*pi*r))*1*k*(lambda/2).*exp(-j*k*r).*sin(theta)];
% efield=abs(efield);
% efield=efield'.*ones(1,length(phi));
% efield=efield';