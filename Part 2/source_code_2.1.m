%WR90
a=22.86*10^(-3);b=10.16*10^(-3);
scdistance=10^(-3)*[48.4 73.0 97.5 122.2];
%Material 1
d1=1.5*10^(-3);
mat1a=10^(-3)*[46.9 71.5 96.0 120.7];
mat1b=10^(-3)*[55.6 80.2 104.7 129.4];
ans1=er(a,scdistance,mat1a,mat1b)
ans1onlya=eronlya(a,scdistance,mat1a-0.0000555,d1)
%Material 2
d2=1.517*10^(-3);
mat2a=10^(-3)*[46.7 71.3 95.8 120.5];
mat2b=10^(-3)*[52.6 77.0 101.7 126.3];
ans2=er(a,scdistance,mat2a,mat2b)
ans2onlya=eronlya(a,scdistance,mat2a+0.0000415,d2)

function er=er(a,d0,d1,d2)
%d0 matrix of short circuit plain
%d1 matrix of first measurement
%d2 matrix of second measurement
c=3*10^8;
fc=c/(2*a);
lg=[0 0 0];
for i=1:3
    lg(i)=d0(i+1)-d0(i);
end
lgmean=2*mean(lg);
f=sqrt(fc^2+(c/lgmean)^2);
%lg=physconst('LightSpeed')/(sqrt(f^2-fc^2));
k=(fc/f)^2;
for i=1:3
    da(i)=d1(i+1)-d0(i);
end
for i=1:3
    db(i)=d2(i)-d0(i);
end
damin=mean(da);
dbmin=mean(db);

er=k-((1-k)/(tan(2*pi*damin/lgmean)*tan(2*pi*dbmin/lgmean)));
end

function eronlya=eronlya(a,d0,d1,d)
%d0 matrix of short circuit plain
%d1 matrix of first measurement
%d=length of dielectric
c=3*10^8;
fc=c/(2*a);
lg=[0 0 0];
for i=1:3
    lg(i)=d0(i+1)-d0(i);
end
lgmean=2*mean(lg);
f=sqrt(fc^2+(c/lgmean)^2);
k=(fc/f)^2;
for i=1:3
    da(i)=d1(i+1)-d0(i);
end
damin=mean(da);
zo=120*pi/sqrt(1-k);
xa=-zo*tan(2*pi*damin/lgmean);
resb=c*xa/(2*pi*120*pi*f*d);
xzero=(2*pi*f*d/c)*sqrt(2-k^2);
x = fzero(@(x) (tan(x)/x-resb) ,xzero);
eronlya=k+(c*x/(2*pi*f*d))^2;
end
