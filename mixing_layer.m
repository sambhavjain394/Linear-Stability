clear all; close all; clc
% Solving for characteristics of mixing layer

% N=input('Enter the value of N: ');
% alpha = input('Enter the value of alpha: ');
% beta=input('Enter the value of beta: ');

N=60;
alpha=0.5;
beta=0;
L_inf=10;
l=1;
[y,DM]=chebdif(N,2);
s=(l/L_inf)^2;

f=@(y) -l*y./sqrt(1+s-y.^2);
f_prime=@(y) -l*(1+s)./(1+s-y.^2).^1.5;
f_dprime=@(y) -3*l*(1+s)*y./(1+s-y.^2).^2.5;
y_mapped=f(y);


% U is the base flow u=0.5*(1+tanh(y))
U = diag(0.5*(1+tanh(y_mapped)));
U_dprime=diag(-1*sech(y_mapped).^2.*tanh(y_mapped));

y_prime=f_prime(y);
y_dprime=f_dprime(y);

y_base_prime=1./y_prime;
y_base_dprime=-y_dprime./(y_prime).^3;

y_dash=diag(y_base_prime);
y_double_dash=diag(y_base_dprime);

D_1=y_dash*DM(:,:,1);
D_2=y_double_dash*DM(:,:,1)+diag(y_base_prime.^2)*DM(:,:,2);
% Rayleigh Equation
k=sqrt(alpha^2+beta^2);
L_os=1i*alpha*U*(k^2*eye(N)-D_2)+1i*alpha*U_dprime;
A=L_os;
B=1i*(k^2*eye(N)-D_2);

%% Applying Boundary Conditions

I=eye(N);
A(1,:)=I(1,:);
A(N,:)=I(N,:);
B([1,N],:)=zeros(2,N);

%% Solving for eigen values 

[V,lambda]=eig(A,B,'qz');
j=1;
for i=1:N
    if ~isinf(imag(lambda(i,i))) && abs(imag(lambda(i,i)))<5 && ~isinf(real(lambda(i,i)))&& real(lambda(i,i))<5
        egn_val(1,j)=lambda(i,i);
        egn_vect(:,j)=V(:,i);
        j=j+1;
    end
end

[im,seq]= sort(imag(egn_val),'ascend');
egn_vect(:,:) = egn_vect(:,seq);
egn_val(:,:)=egn_val(1,seq);

%% Plotting results
figure(1)
plot(real(egn_val),imag(egn_val),'rsq','linewidth',2); 
xlim([0 0.5]);
ylim([-0.1 0.1]);
grid on;

figure(2)
plot(abs(egn_vect(1:N,end)),y_mapped,'b');
grid on;
