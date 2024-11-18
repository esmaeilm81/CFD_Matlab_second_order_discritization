clc;
clear;
close all;

h=[0.25,0.125,0.0625,0.03125,0.015625];
N=(1./h)+1;
N=floor(N);
Norm=zeros(1,length(h));
tol=1e-12;
maxit=N;
Val=zeros(length(h),3);
z=1;

for i=1:length(h)
    A=zeros(N(i),N(i));
    B=zeros(N(i),1);
    for j=1:N(i)
        x(j)=(j-1)*h(i);
        S(j)=x(j)*(1-x(j));
    end
    S=S';
    A(1,1)=1; A(N(i),N(i))=1;
    B(1)=2; B(N(i))=5;
    for j=2:N(i)-1
        B(j)=S(j)*h(i)^2;
    end
    k=2;
    for j=2:N(i)-1
        A(k,j-1)=1;
        A(k,j)=-2;
        A(k,j+1)=1;
        k=k+1;
    end
    A=sparse(A);
    phi=gmres(A,B,[],tol,maxit(i));
    for c=1:3
        Val(i,c)=phi(1+c*z);
    end
    z=z*2;
end
for i=1:length(h)-1
    sum=0;
    for j=1:3
        sum=sum+(Val(i+1,j)-Val(i,j))^2;
    end
    Norm(i)=sqrt(sum);
end
coef=polyfit(log10(h(1:4)),log10(Norm(1:4)),1);
fit=polyval(coef,log10(h(1:4)));
fit=10.^fit;
figure();
loglog(h(1:4),Norm(1:4),'o');
hold on;
loglog(h(1:4),fit,'r');
xlabel('h (log scale)');
ylabel('Norm \phi_i (log scale)');
legend('Data Points','Fitted Line','location','northwest');
