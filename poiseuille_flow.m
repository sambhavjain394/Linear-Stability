clear all; close all; clc

% Solving for characteristics of fully developed channel flow
N=input('Enter the value of N: ');
Re=input('Enter the value of reynolds number: ');
alpha = input('Enter the value of alpha: ');
beta=input('Enter the value of beta: ');

[y,DM]=chebdif(N,4);
% U is the base flow u=1-y^2

U = diag(1-y.^2);
U_prime=diag(-2.*y);
U_dprime=-2*eye(N);
% Linearized Navier Stokes Equation
L=(DM(:,:,2)-(alpha^2+beta^2)*eye(N) )/Re - 1i*alpha*U;
A_1=[1i*alpha*eye(N), DM(:,:,1), 1i*beta*eye(N), zeros(N);
    L, -U_prime, zeros(N), -1i*alpha*eye(N);
    zeros(N), L, zeros(N), -DM(:,:,1);
    zeros(N), zeros(N), L, -1i*beta*eye(N);];

B_1=[zeros(N), zeros(N), zeros(N), zeros(N);
    -1i*eye(N), zeros(N), zeros(N), zeros(N);
    zeros(N), -1i*eye(N), zeros(N), zeros(N);
    zeros(N), zeros(N), -1i*eye(N), zeros(N);];
                       
% Orr-Sommerfield Equation
k=sqrt(alpha^2+beta^2);
L_os=1i*alpha*U*(k^2*eye(N)-DM(:,:,2))+1i*alpha*U_dprime+ (k^4*eye(N)+DM(:,:,4)-2*k^2*DM(:,:,2))/Re;

A_2=L_os;
B_2=1i*(k^2*eye(N)-DM(:,:,2));

% Orr-Sommerfield Equations with Squire Equation
L_sq=1i*alpha*U+(k^2*eye(N)-DM(:,:,2))/Re;
A_3=[L_os, zeros(N), zeros(N);
    1i*beta*U_prime, 1i*beta*L_sq, -1i*alpha*L_sq;
    DM(:,:,1), 1i*alpha*eye(N), 1i*beta*eye(N);];

B_3=[1i*(k^2*eye(N)-DM(:,:,2)), zeros(N), zeros(N);
    zeros(N), -beta*eye(N), alpha*eye(N);
    zeros(N), zeros(N), zeros(N);];

%% Applying Boundary Conditions

% In LNSE
I1=eye(4*N);
A_1([1,N,N+1,2*N,2*N+1,3*N],:)=I1([1,N,N+1,2*N,2*N+1,3*N],:);
B_1([1,N,N+1,2*N,2*N+1,3*N],:)=zeros(6,4*N);

B_1(3*N+1,:)=zeros(1,4*N);
A_1(3*N+1,:)=zeros(1,4*N);
Replacement=DM(:,:,2)-((k^2)*eye(N));
A_1([3*N+1,4*N],3*N+1:4*N)=Replacement([1,N],:);

% In Orr-sommerfield
I2=eye(N);
A_2(1,:)=I2(1,:);
A_2(N,:)=I2(N,:);
A_2(2,:)=DM(1,:,1);
A_2(N-1,:)=DM(N,:,1);
B_2([1,2,N-1,N],:)=zeros(4,N);

% In Orr-Sommerfield equations+ Squirre Equation
I3=eye(3*N);
A_3(1,:)=I3(1,:);
A_3(N+1,:)=I3(N+1,:);
A_3(2*N+1,:)=I3(2*N+1,:);
A_3(N,:)=I3(N,:);
A_3(2*N,:)=I3(2*N,:);
A_3(3*N,:)=I3(3*N,:);

A_3(2,:)=zeros(1,3*N);
A_3(N-1,:)=zeros(1,3*N);
A_3(2,1:N)=DM(1,:,1);
A_3(N-1,1:N)=DM(N,:,1);
B_3([1,2,N-1,N,N+1,2*N,2*N+1,3*N],:)=zeros(8,3*N);

%% Solving for eigen values 

[V1,lambda_1]=eig(A_1,B_1,'qz');
[V2,lambda_2]=eig(A_2,B_2,'qz');
[V3,lambda_3]=eig(A_3,B_3,'qz');
j=1; p=1;r=1;
for i=1:4*N
    if ~isinf(imag(lambda_1(i,i))) && abs(imag(lambda_1(i,i)))<5 && ~isinf(real(lambda_1(i,i)))&& real(lambda_1(i,i))<5
        egn_val_1(1,j)=lambda_1(i,i);
        egn_vect_1(:,j)=V1(:,i);
        j=j+1;
    end
    if i<=N && ~isinf(imag(lambda_2(i,i))) && abs(imag(lambda_2(i,i)))<5 && ~isinf(real(lambda_2(i,i)))&& real(lambda_2(i,i))<5
        egn_val_2(1,p)=lambda_2(i,i);
        egn_vect_2(:,p)=V2(:,i);
        p=p+1;
    end
    if i<=3*N && ~isinf(imag(lambda_3(i,i))) && abs(imag(lambda_3(i,i)))<5 && ~isinf(real(lambda_3(i,i)))&& real(lambda_3(i,i))<5
        egn_val_3(1,r)=lambda_3(i,i);
        egn_vect_3(:,r)=V3(:,i);
        r=r+1;
    end
end

[im_1,seq_1]= sort(imag(egn_val_1),'ascend');
egn_vect_1(:,:) = egn_vect_1(:,seq_1);
egn_val_1(:,:)=egn_val_1(1,seq_1);

[im_2,seq_2]= sort(imag(egn_val_2),'ascend');
egn_vect_2(:,:) = egn_vect_2(:,seq_2);
egn_val_2(:,:)=egn_val_2(1,seq_2);

[im_3,seq_3]= sort(imag(egn_val_3),'ascend');
egn_vect_3(:,:) = egn_vect_3(:,seq_3);
egn_val_3(:,:)=egn_val_3(1,seq_3);


%% Plotting results
figure(1)
plot(real(egn_val_1),imag(egn_val_1),'g*'); 
hold on;
plot(real(egn_val_2),imag(egn_val_2),'sq');
plot(real(egn_val_3),imag(egn_val_3),'r+');
xlim([0 1]);
ylim([-1 0.4]);
grid on;

figure(2)
plot(y(:,1),abs(egn_vect_1(1:N,end)),'b');
hold on;
plot(y(:,1),abs(egn_vect_1(N+1:2*N,end)),'r');
plot(y(:,1),abs(egn_vect_1(2*N+1:3*N,end)),'g');
plot(y(:,1),abs(egn_vect_1(3*N+1:4*N,end)),'k');
grid on;
