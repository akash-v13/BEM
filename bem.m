clear
NO=input('Enter the no. of boundary elements per side:')
N=4*NO;
dl=1/NO;
for i=1:NO
    xb(i)=(i-1)*dl;
    yb(i)=0;
    xb(NO+i)=1;
    yb(NO+i)=xb(i);
    xb(2*NO+i)=1-xb(i);
    yb(2*NO+i)=1;
    xb(3*NO+i)=0;
    yb(3*NO+i)=1-xb(i);
end
xb(N+1)=xb(1);
yb(N+1)=yb(1);
for i=1:N
    xm(i)=0.5*(xb(i)+xb(i+1));
    ym(i)=0.5*(yb(i)+yb(i+1));
    lg(i)=sqrt((xb(i+1)-xb(i))^2+(yb(i+1)-yb(i))^2);
    nx(i)=(yb(i+1)-yb(i))/lg(i);
    ny(i)=(-xb(i+1)+xb(i))/lg(i);
end
for i=1:N
    if i<=NO
        BCT(i)=1;
        BCV(i)=0;
    elseif i>NO & i<=2*NO
        BCT(i)=0;
        BCV(i)=cos(pi*ym(i));
    elseif i>2*NO & i<=3*NO
        BCT(i)=1;
        BCV(i)=0;
    else
        BCT(i)=0;
        BCV(i)=0;
    end
end
for m=1:N
    BC(m)=0;
    for k=1:N
        A(k)=lg(k)^2;
        B(k)=2*lg(k)*(-ny(k)*(xb(k)-xm(m))+nx(k)*(yb(k)-ym(m)));
        E(k)=(xb(k)-xm(m))^2+(yb(k)-ym(m))^2;
        D(k)=sqrt(abs(4*A(k)*E(k)-B(k)^2));
        BA(k)=B(k)/A(k);
        EA(k)=E(k)/A(k);
        if 4*A(k)*E(k)-B(k)^2 == 0
            PF1(k)=0.5*lg(k)*(log(lg(k))+(1+0.5*BA(k))*log(abs(1+0.5*BA(k)))-0.5*BA(k)*log(abs(0.5*BA(k)))-1);
            PF2(k)=0;
        else
            PF1(k)=0.25*lg(k)*(2*(log(lg(k))-1)-0.5*BA(k)*log(abs(EA(k)))+(1+0.5*BA(k))*log(abs(1+BA(k)+EA(k)))+(D(k)/A(k))*(atan((2*A(k)+B(k))/D(k))-atan(B(k)/D(k))));
            PF2(k)=lg(k)*(nx(k)*(xb(k)-xm(m))+ny(k)*(yb(k)-ym(m)))/D(k)*(atan((2*A(k)+B(k))/D(k))-atan(B(k)/D(k)));
        end
        
        F1(k)=PF1(k)/pi;
        F2(k)=PF2(k)/pi;
        
        if k==m
            del=1;
        else
            del=0;
        end
        if BCT(k)==0
            AB(m,k)=-F1(k);
            BC(m)=BC(m)+BCV(k)*(-F2(k)+0.5*del);
        else
            AB(m,k)=F2(k)-0.5*del;
            BC(m)=BC(m)+BCV(k)*F1(k);
            
        end
        
    end
end
NAC = max(size(AB));
% Perform Gaussian Elimination
 for j=2:NAC,
     for i=j:NAC,
        m = AB(i,j-1)/AB(j-1,j-1);
        AB(i,:) = AB(i,:) - AB(j-1,:)*m;
        BC(i) = BC(i) - m*BC(j-1);
     end
 end
% Perform back substitution
 Z = zeros(NAC,1);
 Z(NAC) = BC(NAC)/AB(NAC,NAC);
 for j=NAC-1:-1:1,
   Z(j) = (BC(j)-AB(j,j+1:NAC)*Z(j+1:NAC))/AB(j,j);
 end
for m=1:N
    if BCT(m)==0
        phi(m)=BCV(m);
        dphi(m)=Z(m);
    else
        phi(m)=Z(m);
        dphi(m)=BCV(m);
    end
end
xi=input('Enter the point xi:')
eta=input('Enter the point eta:')
sum=0;
for i=1:N
        A(i)=lg(i)^2;
        B(i)=2*lg(i)*(-ny(i)*(xb(i)-xi)+nx(i)*(yb(i)-eta));
        E(i)=(xb(i)-xi)^2+(yb(i)-eta)^2;
        D(i)=sqrt(abs(4*A(i)*E(i)-B(i)^2));
        BA(i)=B(i)/A(i);
        EA(i)=E(i)/A(i);
        if 4*A(i)*E(i)-B(i)^2 == 0
            PF1(i)=0.5*lg(i)*(log(lg(i))+(1+0.5*BA(i))*log(abs(1+0.5*BA(i)))-0.5*BA(i)*log(abs(0.5*BA(i)))-1);
            PF2(i)=0;
        else
            PF1(i)=0.25*lg(i)*(2*(log(lg(i))-1)-0.5*BA(i)*log(abs(EA(i)))+(1+0.5*BA(i))*log(abs(1+BA(i)+EA(i)))+(D(i)/A(i))*(atan((2*A(i)+B(i))/D(i))-atan(B(i)/D(i))));
            PF2(i)=lg(i)*(nx(i)*(xb(i)-xi)+ny(i)*(yb(i)-eta))/D(i)*(atan((2*A(i)+B(i))/D(i))-atan(B(i)/D(i)));
        end
sum=sum+phi(i)*PF2(i)-dphi(i)*PF1(i);
end            
pint=sum/pi