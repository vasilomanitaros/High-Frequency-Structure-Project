N=201;                  %ορισμός τιμών συχνοτήτων
f0=10^9;                %ορισμός συχνότητας αναφοράς για τα μήκη
f=0:(4*f0/N):(4*f0);    %ορισμός διανύσματος συχνοτήτων
zfreqfortio=100+1./(j*2*pi*2*10^(-12).*f);  %διάνυσμα χωρητικότητας φορτίου με τη συχνότητα
zfreqpyknwtis=1./(j*2*pi*2.7*10^(-12).*f);  %διάνυσμα χωρητικότητας πυκνωτή στην είσοδο με τη συχνότητα
zin1=zin(50,zfreqfortio,0.2,f,f0);  %υπολογισμός αντίστασης εισόδου από το φορτίο μέχρι την είσοδο
zin2=zin(50,0,0.13,f,f0);
zina=parallila(zin1,zin2);
zinb=zin(50,zina,0.1,f,f0);
zintotal=parallila(zinb,zfreqpyknwtis);
reflgamma=(zintotal-50)./(zintotal+50); %υπολογισμός συντελεστή ανάκλασης ως μιγαδικός
abgamma=abs(reflgamma); %μέτρο συντελεστή ανάκλασης
loggamma=20*log10(abgamma); %μέτρο σε dB
figure(1);
plot(f,abgamma);
ans11=interp1(f,abgamma,10^9);%εύρεση μέτρου συντελεστή ανάκλασης για το ερώτημα 1.1α
ans11b=interp1(f,abgamma,1.5*10^9);%εύρεση μέτρου συντελεστή ανάκλασης για το ερώτημα 1.1b
figure(2);
plot(f,loggamma);

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
