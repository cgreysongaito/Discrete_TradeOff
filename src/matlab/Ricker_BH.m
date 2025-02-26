%orbit diagram Ricker & constant

a = 2;
b = 1;
pbar = 0.4; %0.8;
alpha = 0.01; %608102318/8513732705;
K = 1;
beta = 0.3;

tau = 2;

%analytical check q<1?
q=(a-b*exp(-K*(tau+1)))*pbar^(tau+1);
disp(strcat('for solutions to be bounded I need ', num2str(q), ' to be less than 1'));

taur = 1/K*log(b*(K+log(1/pbar)/(a*log(1/pbar))))-1;
%taur<0 -> peak of q appeared before 0
% q(tau=0)>1-exp(-alpha)
qr=(a-b*exp(-K*(taur+1)))*pbar^(taur+1);
%disp('is qr less than 1-exp(-alpha)?');
%disp(qr<1-exp(-alpha));
disp(qr-(1-exp(-alpha)));

%Ricker-constant
%Nv = [10*rand(1, tau+1)];
Nv = [3.7957    7.0405    1.9182];

for i=1:50000
    Nnew = exp(-alpha-beta*Nv(end))*Nv(end)+(a-b*exp(-K*(tau+1)))*pbar^(tau+1)*Nv(end-tau);
    Nv = [Nv, Nnew];
 

end


plot(Nv,'*');
xlabel('time');
ylabel('sol-Ricker constant')



