function  [lambda0,n]=  Intersect_cylinder(a,v,b,w,R)
% Berekent de kleinste positieve waarde lambda0 zodat de lijn a+lambda*v voor lambda=lambda0>0 de cilinder met
% hartlijn b+mu*w en straal R doorsnijdt, n is dan de normaal, gericht naar buiten, op de cilinder in het snijpunt
% Als er geen snijpunt bestaat voor een lambd0>0, dan wordt een groot positief getal geretourneerd
% en voor n de nulvector.
% Note: lambda0>0, hence the intersection point is NOT the start point a of the line
% eerst v en w normaliseren
vn=v/norm(v);
wn=w/norm(w);
A=1-(vn'*wn)^2;
B=2*(a-b)'*(vn-(vn'*wn)*wn);
C=(a-b)'*(a-b)-((a-b)'*wn)^2-R^2;

if A==0
     % omdat vn en wn beiden genormaliseerd zijn , kan dit alleen als vn=wn of vn=-wn
     % dan zijn lijn en hartlijn van de cilinder evenwijdig
     lambda0=10^15; 
     n=[0;0;0];
else
     % nu echte kwadratische vergelijking
     Det=B^2-4*A*C;
     if Det<0 % geen reële oplossingen
         lambda0=10^15;
     else %Det >=0;
         lambda1=(-B+sqrt(Det))/(2*A);
         lambda2=(-B-sqrt(Det))/(2*A);
         % omdat A>0 is altijd lambda1>=lambda2
         if lambda2>1e-10   %lambda2 is de kleinste positive oplossing
              lambda0=lambda2;
         elseif lambda1>1e-10 
              % lambda2<=0 en lambda1>0, dus lambda1 is de kleinste positive oplossing
              lambda0=lambda1;
         else %% lambda1<=0 en lambda2<=0, geen positieve oplossingen
              lambda0=10^15;
         end;
    end;
end;
n=[0;0;0];
if lambda0<10^15
    % er was echt snijpunt
    x= a+lambda0*vn;      % x=snijpunt lijn met cilinder
    y= b+((x-b)'*wn)*wn;  % y=bijbehorend punt op hartlijn van cilinder  
    n=(x-y)/(norm(x-y));   % n=normaal op cilinder in snijpunt, gericht naar buiten
    lambda0=lambda0/norm(v); 
    %  het snijpunt is nu a+lambda0*v, ook als v niet genormaliseerd was 
end;
