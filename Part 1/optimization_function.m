function optgamma=optgamma(p)
z0=50;
normf = (0.01:0.01:2);
zload=10+j*15;
z1=zin(z0,zload,p(1),normf);
zo1=zin(z0,inf,p(4),normf);
za=parallila(zo1,z1);
z2=zin(z0,za,p(2),normf);
zo2=zin(z0,inf,p(5),normf);
zb=parallila(zo2,z2);
z3=zin(z0,zb,p(3),normf);
zo3=zin(z0,inf,p(6),normf);
ztot=parallila(zo3,z3);
reflgamma=(ztot-z0)./(ztot+z0);
abgamma=abs(reflgamma);
loggamma=20*log10(abgamma);
%figure(1)
 %plot(normf,abgamma)
plot(normf,loggamma)
lbest=0;
for i=1:length(normf)
    l=0;
    k=0;
    while (loggamma(i+k)<=-10 && i+k<length(normf))
        l=l+1;
        k=k+1;
    end
    if (l>lbest)
        lbest=l;
    end
end
x=2*lbest/200;
optgamma=2-x;
function zin= zin(z0,zl,l,normf)
    if zl==Inf
        zin=-j*z0.*cot(2*pi*l.*normf)
    else
        zin=z0.*(zl+j*z0.*tan(2*pi*l.*normf))./(z0+j*zl.*tan(2*pi*l.*normf));
    end
end

    function parallila= parallila(z1,z2)
        y1=1./z1;
        y2=1./z2;
        ypar=y1+y2;
        parallila=1./ypar;
    end
end