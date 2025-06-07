t1_values=[1,1.20185];t2=1;gamma=4./3;
figure;hold on;
J1=[];J2=[];
for k=1:2
    t1=t1_values(k);
    a=t2.*(t1+gamma./2);
    c=t2.*(t1-gamma./2);
    E_values=-3:0.001:3;
    for j=1:length(E_values)
        E=E_values(j);
        b=(t1.^2-gamma.^2./4)+t2.^2-E.^2;
        beta=roots([a,b,c]);
        J1=[J1;abs(beta)];
        J2=[J2;abs(E)*ones(2,1)];
        % 同样为了画图时对齐
    end

    scatter(J2,J1,1,"black",'filled');
end