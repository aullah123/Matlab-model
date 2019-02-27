clear all
close all
clc


n =10;
k = 1;
h=2;

xo=5;
yo=5;

x=0:10;
y=0:10;

for i=1:n+1
    for j=1:n+1
        alpha(i,j) = 0.01;
        beta(i,j) = 0.3;
        gamma(i,j) = 0.5;
    end
end

delta1=0.1;
delta2=0.2;
delta3=0.3;

for i=1:n+1
    for j=1:n+1
        S0(i,j)=2*exp(-((i-xo).^2 + (j-yo).^2  )); 
        A0(i,j)=1*exp(-((i-xo).^2 + (j-yo).^2  )); %initial addicted people
        R0(i,j)=0;
    end
end
S=S0;
Snew=S0;
A=A0;
Anew=A0;
R=R0;
Rnew=R0;
for m=1:20
    for i=2:n
        for j=2:n
            Snew(i,j)=S(i,j)+k*(-alpha(i,j)*S(i,j)*A(i,j) ...
                +delta1*((S(i-1,j)-2*S(i,j)+S(i+1,j))+S(i,j-1)-2*S(i,j)+S(i,j+1))/h^2)
            Anew(i,j)=A(i,j)+k*(alpha(i,j)*S(i,j)*A(i,j) ...
                +delta2*((A(i-1,j)-2*A(i,j)+A(i+1,j))+A(i,j-1)-2*A(i,j)+A(i,j+1))/h^2 ...
                +gamma(i,j)*R(i,j)*A(i,j));
            Rnew(i,j)=R(i,j)+k*(beta(i,j)*A(i,j)+delta3*(((R(i-1,j)-2*R(i,j)...
                +R(i+1,j))/h^2)+((R(i,j-1)-2*R(i,j)+R(i,j+1))/(h^2)))-gamma(i,j)*R(i,j)*A(i,j));
        end
    end
    for j=2:n
        Snew(1,j)=0;      %left boundary
        Snew(n+1,j)=0;    %right boundary
        Anew(1,j)=0;      %left boundary
        Anew(n+1,j)=0;    %right boundary
        Rnew(1,j)=0;      %left boundary
        Rnew(n+1,j)=0;    %right boundary
    end
    for i=2:n
        Snew(i,1)=0;      %bottom boundary
        Snew(i,n+1)=0;    %top boundary 
        Anew(i,1)=0;      %bottom boundary
        Anew(i,n+1)=0; 
        Rnew(i,1)=0;      %bottom boundary
        Rnew(i,n+1)=0; 
    end
    subplot(1,3,1)
    surf(x,y,Snew);
    title(['time: ',num2str(m*k),' seconds'])
    shading interp;
    zlim([0 2])
    drawnow
    subplot(1,3,2)
    surf(x,y,Anew);
    title(['time: ',num2str(m*k),' seconds'])
    shading interp;
    zlim([0 2])
    drawnow
    subplot(1,3,3)
    surf(x,y,Rnew);
    title(['time: ',num2str(m*k),' seconds'])
    shading interp;
    zlim([0 2])
    drawnow
    F(m) = getframe;
    S=Snew;
    A=Anew;
    R=Rnew;
end
movie2avi(F,'test.avi')