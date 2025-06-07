
t2=1;
gamma=4/3;
t1_values=linspace(-2, 2, 1000); 
critical_point=sqrt(t2^2+(gamma/2)^2); 


W_values=zeros(size(t1_values));
for i=1:length(t1_values)
 t1=t1_values(i);
 t1_tilde = sqrt((t1-gamma/2)*(t1+gamma/2));
 if abs(t1_tilde)<t2
        W=1;
    else
        W=0;
    end
    W_values(i)=W;
end

figure;
plot(t1_values, W_values, 'b-', 'LineWidth', 2);
hold on;

grid on;