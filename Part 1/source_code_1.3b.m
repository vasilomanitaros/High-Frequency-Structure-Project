N=201;          %oρισμός τιμών συχνοτήτων
f0=10^9;        %ορισμός συχνότητας αναφοράς για τα μήκη
f = 0:(3*f0/N):(3*f0);  %ορισμός διανύσματος συχνοτήτων
zload=50;   %επίλυση κυκλώματος
z1=zin(98.45,inf,1/8,f,f0);
z1r=parallila(z1,zload);
z2=zin(101.6,z1r,1/8,f,f0);
z3=zin(43.6,inf,1/8,f,f0);
za=parallila(z2,z3);
z4=zin(101.6,za,1/8,f,f0);
z5=zin(98.45,inf,1/8,f,f0);
zb=parallila(z4,z5);
reflgamma=(zb-50)./(zb+50); %υπολογισμός συντελεστή ανάκλασης
abgamma=abs(reflgamma);
loggamma=20*log10(abgamma);
swr=(1+abgamma)./(1-abgamma);   %υπολογισμός SWR
for k=1:201 %εφαρμογή περιορισμού για το SWR
    if swr(k)>10;
        swr(k)=10;
    end
end
figure(2)
plot(f,swr)
for k=1:201 %εφαρμογή περιορισμού για τον συντελεστή ανάκλασης
    if loggamma(k)<-60;
        loggamma(k)=-60;
    end
end
figure(1)
plot(f,loggamma)
loggamma1=loggamma(66:101); %υπολογισμός συχνοτήτων αποκοπής (|Γ|)db> -10dB
loggamma2=loggamma(169:202);
f1=f(66:101);
f2=f(169:202);
flow=interp1(loggamma1,f1,-10);
fhigh=interp1(loggamma2,f2,-10);

function zin= zin(z0,zl,l,f,f0)
    if zl==Inf
        zin=-j*z0.*cot((2*pi*l/f0).*f)
    else
        zin=z0.*(zl+j*z0.*tan((2*pi*l/f0).*f))./(z0+j*zl.*tan((2*pi*l/f0).*f));
    end
end

function parallila=parallila(z1,z2)
    y1=1./z1;
    y2=1./z2;
    ypar=y1+y2;
    parallila=1./ypar;
end