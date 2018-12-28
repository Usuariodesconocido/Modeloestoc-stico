% Este programa resuelve num?ricamente el modelo de crecimiento de un
% sector con utilidad logaritmica y funci?n de producci?n Cobb-Douglas
% calibrado en la clase



%TARGETS DE CALIBRACI?N
ke=2.9;
ie=0.1671;
he=0.25;
ce=1-ie;
%VALORES DE LOS PAR?METROS
% Tasa crecimiento de la oferta de trabajo 
n=0.0173;

% Tasa crecimiento de la productividad 
ga=0.0139;

% Tecnolog?a 
d=(ie/ke)+1-((1+ga)*(1+n));
a=0.3;
th=(ke^(-a))*(he^(a-1));

% Preferencias
b=((1+ga)*(1+n))/((a*th*(ke^(a-1))*(he^(1-a)))+1-d);
sg=(1-he)*((1-a)*th*(ke^a)*(he^(-a)))/ce;


% Proceso estoc?stico para TFP

rho=0.95;
mu=(1-rho)*th;
sig=2*mu;
q=85;
[ths,Pi]=mytauchen(mu,rho,sig,q);

% Malla de valores para k
Kmax=4;
Kmin=0.1;
p=200;
eta=(Kmax-Kmin)/(p-1);
K(1)=Kmin;
for i=2:p
    K(i)=Kmin+((i-1)*eta);
end

% Malla de valores para h
Hmax=0.7;
Hmin=0.06;
ph=200;
etah=(Hmax-Hmin)/(ph-1);
H(1)=Hmin;
for i=2:ph
    H(i)=Hmin+((i-1)*etah);
end

% La matriz M

for i=1:p
    for m=1:q
        
    for h=1:ph
        for j=1:p
        if (ths(m)*(K(i)^a)*(H(h)^(1-a)))+((1-d)*K(i))-((1+n)*(1+ga)*K(j))>0 
            M(p*(m-1)+i,p*(h-1)+j)=log((ths(m)*(K(i)^a)*(H(h)^(1-a)))+((1-d)*K(i))-((1+n)*(1+ga)*K(j)))+(sg*log(1-H(h)));
        else
            M(p*(m-1)+i,p*(h-1)+j)=-10000;
        end
        end
    end
    
    end
end

% La funci?n de valor inicial V0 

V0=zeros(p*q,1);
v=V0;
V=ones(p*q,1);

% Iterando la funci?n de valor
while norm(V-v)>10^(-6)    
    v=V;
    for i=1:q
        for s=1:p
            W(i,s)=0;
            for j=1:q            
                W(i,s)=W(i,s)+(Pi(i,j)*V(((j-1)*p)+s));
            end
        end
    end
    
for i=1:p
    for m=1:q
        
    for h=1:ph
        for j=1:p
                BEL(p*(m-1)+i,p*(h-1)+j)=M(p*(m-1)+i,p*(h-1)+j)+(b*W(m,j));
        end
    end
    
    end
end
    [V g]=max(BEL,[],2);
end
        
        
% % La funci?n de valor v(k)
for i=1:q
    for j=1:p
        VV(i,j)=V(p*(i-1)+j);
    end
end
[X,Y]=meshgrid(K,ths);
%mesh(X,Y,VV)
surf(X,Y,VV,'FaceColor','red','EdgeColor','none')
camlight left; 
lighting phong
pause

 % La funci?n de pol?tica para k
 for i=1:q
     for j=1:p
         h=ceil(g(((i-1)*p)+j)/p);
         s=g(((i-1)*p)+j)-(p*(h-1));
         K1(i,j)=K(s);
         H1(i,j)=H(h);
     end
 end

%mesh(X,Y,K1)
surf(X,Y,K1,'FaceColor','red','EdgeColor','none')
camlight left; 
lighting phong
pause        
        
% La funci?n del control h        
        
surf(X,Y,H1,'FaceColor','red','EdgeColor','none')
camlight left; 
lighting phong
pause     
        




        
        