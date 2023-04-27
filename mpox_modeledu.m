function [S] = mpox_modeledu(par, time);

alpha=par(1);
beta=par(2);
beta_e=par(3);
sigma=par(4);
gamma=par(5);
S0=par(6);
S_e0=par(7);
E0=par(8);
E_e0=par(9);
I0=par(10);
I_e0=par(11);
R0=par(12);

f = @(t,x) [-beta*x(1)*x(5)-alpha*x(1); 
    alpha*x(1)-beta_e*x(2)*x(5);
    beta*x(1)*x(5)-alpha*x(3)-sigma*x(3);
    alpha*x(3)+beta_e*x(5)*x(2)-sigma*x(4);
    sigma*x(3)-alpha*x(5)-gamma*x(5);
    alpha*x(5)+sigma*x(4)-gamma*x(6);
    gamma*x(5)+gamma*x(6)];

[~, Y] = ode45(f, 0:1:time, [S0,S_e0,E0,E_e0,I0,I_e0,R0]);

out = zeros(7,length(0:1:time));
out(1,:) = (Y(:,1))';
out(2,:) = (Y(:,2))';
out(3,:) = (Y(:,3))';
out(4,:) = (Y(:,4))';
out(5,:) = (Y(:,5))';
out(6,:) = (Y(:,6))';
out(7,:) = (Y(:,7))';

S = round(out(1,1),9)-round(out(1,150),9)+round(out(2,1),9)-round(out(2,150),9);

plot(0:1:time, Y(:,5) + Y(:,6), 'LineWidth', 1.5);
xlabel('Days');
ylabel('Population Infected');
ylim([0, 3000]);
legend('None','0.1% per day','0.2% per day','0.3% per day', '0.4% per day');
title("Total Infected by Educational Coverage");


end

