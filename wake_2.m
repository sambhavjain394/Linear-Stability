clear all; close all; clc

Re=100;
N=60;
alpha=1.4;
beta=0;
L_inf=9;
k=sqrt(alpha^2+beta^2);
[y,DM]=chebdif(N,4);

f=@(y) (y+1)*L_inf/2;
f_prime=@(y) y*0+L_inf/2;
y_mapped=f(y);
y_prime=f_prime(y);
y_dprime=zeros(N,1);
y_tprime=zeros(N,1);
y_qprime=zeros(N,1);

%Base flow
U=diag(1-exp(-(y_mapped.^2)*log(2) ));
U_prime=diag(2*log(2)*y.*exp(-(y_mapped.^2)*log(2)) );
U_dprime=diag(2*log(2)*exp(-(y_mapped.^2)*log(2)).*(1-2*y_mapped.^2*log(2)) );

y_base_prime=1./y_prime;
y_base_dprime=-y_dprime./(y_prime.^3);
y_base_tprime=(3*y_dprime.^2-y_tprime.*y_prime)./(y_prime.^5);
y_base_qprime=(10*y_prime.*y_dprime.*y_tprime-15*y_dprime.^3-y_qprime.*y_prime.^2)./(y_prime.^7);

y_dash=diag(y_base_prime);
y_double_dash=diag(y_base_dprime);
y_triple_dash=diag(y_base_tprime);
y_quad_dash=diag(y_base_qprime);
D1=y_dash*DM(:,:,1);
D2=y_double_dash*DM(:,:,1)+y_dash*y_dash*DM(:,:,2);
D3=y_triple_dash*DM(:,:,1)+3*y_double_dash*y_dash*DM(:,:,2)+y_dash*y_dash*y_dash*DM(:,:,3);
D4=y_quad_dash*DM(:,:,1)+4*y_quad_dash*y_dash*DM(:,:,2)+3*y_double_dash*y_double_dash*DM(:,:,2)+6*y_double_dash*y_dash*y_dash*DM(:,:,3)+y_dash*y_dash*y_dash*y_dash*DM(:,:,4);

% Orr Sommerfield Equation
A=1i*alpha*U*(k^2*eye(N)-D2)+1i*alpha*U_dprime+(k^2*eye(N)-D2)*(k^2*eye(N)-D2)/Re;
B=1i*(k^2*eye(N)-D2);

%% Applying Boundary Conditions
I=eye(N);
A(1,:)=I(1,:);
A(2,:)=D1(1,:);
B([1,2,N-1,N],:)=zeros(4,N);
A_s=A;
A_v=A;
% For sinuous mode
A_s(N,:)=D1(N,:);
A_s(N-1,:)=D3(N,:);
% For varicose mode
A_v(N,:)=I(N,:);
A_v(N-1,:)=D2(N,:);

%% Solving for eigen values 
    
[V_s,lambda_s]= eig( A_s,B,'qz');
[V_v,lambda_v]= eig( A_v,B,'qz');

j=1; l=1;
for i=1:N
    if ~isinf(imag(lambda_s(i,i))) && ~isinf(real(lambda_s(i,i)))
        egn_val_s(1,j)=lambda_s(i,i);
        egn_vect_s(:,j)=V_s(:,i);
        j=j+1;
    end
    if ~isinf(imag(lambda_v(i,i)))&& ~isinf(real(lambda_v(i,i)))
        egn_val_v(1,l)=lambda_v(i,i);
        egn_vect_v(:,l)=V_v(:,i);
        l=l+1;
    end    
end

[im_s,seq_s]= sort(imag(egn_val_s),'ascend');
egn_vect_s(:,:) = egn_vect_s(:,seq_s);
egn_val_s(:,:)=egn_val_s(1,seq_s);

[im_v,seq_v]= sort(imag(egn_val_v),'ascend');
egn_vect_v(:,:) = egn_vect_v(:,seq_v);
egn_val_v(:,:)=egn_val_v(1,seq_v);

plot(y_mapped,abs(egn_vect_v(1:N,end)),'b*');
hold on;
plot(y_mapped,abs(egn_vect_s(1:N,end)),'r+');
legend('Varicose Mode', 'Sinuous Mode');
grid on;
xlim([0 8]); ylim([0 1]);
