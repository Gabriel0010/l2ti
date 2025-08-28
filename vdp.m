function Q=vdp(ImRef,ImDeg)
%ImRef et ImDeg entre 0 et 1

[cortex,base]=Dalycortex;
Lumimag=DalyLum(ImRef);
imagettes=FiltrageCortex(Lumimag,cortex,base,256);
bk=mean(mean(imagettes(1,:,:)));   % la luminance moyenne de l image de base

beta=2; k1=exp(-7*log(6)/3); k2=exp(10*log(6)/3); s=0.7; b=4;

for i=1:31,
    for x=1:256,
        for y=1:256,
           
           if imagettes(i,x,y)==0
                tkl(i,x,y)=1;
            else
                 a=abs(double(imagettes(i,x,y)));
                 a=k2*a;
                 a=exp(s*log(a));
                 a=k1*a;
                 a=exp(b*log(a)); 
                 a=a+1;
                 a=exp(log(a)/b); 
                 tkl(i,x,y)=a;
           end
       end
   end
end

Lumimag=DalyLum(ImDeg);
imagettes2=FiltrageCortex(Lumimag,cortex,base,256);

ckl=abs((imagettes- imagettes2)/bk);

for i=1:31,
    for x=1:256,
        for y=1:256,
          
           if imagettes2(i,x,y)==0
                tkl(i,x,y)=1;
            else
                 a=abs(double(imagettes2(i,x,y)));
                 a=k2*a;
                 a=exp(s*log(a));
                 a=k1*a;
                 a=exp(b*log(a)); 
                 a=a+1;
                 a=exp(log(a)/b);     
           end
           
           tkem=a;   % TKem est la sous images d augmentation du seuil
           if tkl(i,x,y)<a tkem=tkl(i,x,y); end
           if tkem==0  
                   if ckl(i,x,y)==0  
                             Pkl(i,x,y)=0;  
                   else
                             Pkl(i,x,y)=1;
                   end
           else 
               if ckl(i,x,y)==0 
                       Pkl(i,x,y)=0;
               else
                       a=ckl(i,x,y)/tkem;  
                       a=exp(beta*log(a)); 
                       Pkl(i,x,y)=1-exp(-1*a); 
               end
          end
       
       
     end   
   end
end




for k=1:9,
    Seuil(k)=0.1*k;
    n(k)=0;
end

      for x=1:256
            for y=1:256,
               a=1;
                for i=1:31,
                   a=a*(1-Pkl(i,x,y));
               end
               VDP(x,y)=1-a;
               for k=1:9,
                   if (VDP(x,y)>Seuil(k)) | (VDP(x,y)==Seuil(k))  
                   n(k)=n(k)+1;  
                   end
              end
              
           end
       end

Q=n;






function [cortex,base]=Dalycortex;
% function [cortex,base]=Dalycortex;
% cortex(k,l,x,y) k: dom l: fan  Les 30 filtres cortex
% base: le filtre bande de base.
% n: dimension de l'image nombre de pixels
%
%
%
% voir aussi:
%           FiltrageCortex    Decompositioncortex    Dalycsf    DalyLuminance   Daly   ListeDaly   
n=256;
for x=1:n,
    for y=1:n,
   x1=x-((n+1)/2);
   y1=y-((n+1)/2);
         % calcul de la fréquence cycles/pixel
        f(x,y)=sqrt(((x1)*(x1))+((y1)*(y1)))/((n-1)*sqrt(2)); 
        
        if (x1==0)
            if y1>0 
                angle=90;
            else
                angle=270;
            end
        else
        angle=180*atan(y1/x1)/pi;
        end
       teta(x,y)=angle; 
        
    end
end

fk(1)=double(1); 
tw(1)=double(2*fk(1)/3);
for i=2:5,
    fk(i)=fk(i-1)/2;
    tw(i)=2*fk(i)/3;
end

% filtre Bande de Base: BaseBand
fkk=fk(5)/2;
tww=2*fkk/3;
segma=(fkk+(tww/2))/3;
b=fkk+(tww/2);
for x=1:n,
    for y=1:n,
     if f(x,y)< b 
         base(x,y)=double(exp(-(f(x,y)*f(x,y))/(2*segma*segma)));
     else 
         base(x,y)=0;
     end
   end
end
% filtres Dom: 5 filtres
a=fk(5)-(tw(5)/2);
b=fk(5)+(tw(5)/2);
for x=1:n,
    for y=1:n,
        if f(x,y) < a 
            fan(5,x,y)=1; 
        end
        if f(x,y) > b 
            fan(5,x,y)=0; 
        end
        if (f(x,y)>a)&(f(x,y)<b) 
            fan(5,x,y)=(1+cos(pi*((f(x,y)-a)/tw(5))))/2;
        end
        dom(5,x,y)=fan(5,x,y)-base(x,y);
    end
end
% 
for i=1:4,
a=fk(5-i)-(tw(5-i)/2);
b=fk(5-i)+(tw(5-i)/2);

for x=1:n,
    for y=1:n,
        if f(x,y) < a 
            fan(5-i,x,y)=1; end
        if f(x,y) > b 
            fan(5-i,x,y)=0; end
        if (f(x,y)>a)&(f(x,y)<b) 
            fan(5-i,x,y)=(1+cos(pi*((f(x,y)-a)/tw(5-i))))/2;
        end
     dom(5-i,x,y)=fan(5-i,x,y)-fan(6-i,x,y);
    end
end
end

% filtres FAN
delta=30;
for l=1:6,
    tetac=double(((l-1)*delta)-90);
    for x=1:n,
        for y=1:n,
            a=abs(teta(x,y)-tetac);
            if a > delta
                fan(l,x,y)=0;
            else
                fan(l,x,y)=double((1+cos(pi*(a/delta)))/2);
            end
        end
    end
end 
 % Filtres Cortex

 for k=1:5,
     for l=1:6,
        for x=1:n,
            for y=1:n,   
                cortex(k,l,x,y)=fan(l,x,y)*dom(k,x,y);    
            end
        end
    end
end

    


function [Lumimag]=DalyLum(im);
% function [Lumimag]=Dalyluminance(nom);
% nom: fichier image format BMP.
% r/rmax=l/(l+(12.6*l)^0.63)
% l=0.3081*ng+0.1646
%
%
%
% voir aussi:
%           FiltrageCortex    Dalycortex    Dalycsf    Decompositioncortex   Daly   ListeDaly   
im=256*im;
n=size(im);
dimx=n(1);
dimy=n(2);
rmax=0;
for x=1:dimx,
    for y=1:dimy,
        l=0.3081*double(im(x,y))+0.1646;
        r(x,y)=l/(l+exp(0.63*log(12.6*l)));
        if r(x,y)>rmax rmax=r(x,y); end
    end
end
Lumimag=r/rmax;

function [csf]=Dalycsf(n);
% function [csf]=Dalycsf(n);
% n: nombre de pixels
%
%
% voir aussi:
%           FiltrageCortex    Dalycortex    Decompositioncortex    DalyLuminance   Daly   ListeDaly   
p=250;
wm=0.076;   % la taille de l image en métres: 256x256 ---> 7.6 cm.
d=6*wm;     % la distance de visualisation égale six fois la taille de l image.
w=360*atan(1/12)/pi;   % la taille de l image en degré visuel.
la=100;      % la luminance d adaptation (cd/m^2)
ec=0;
al=0.801*exp(-0.2*log(1+(0.7/la)));
bl=0.3*exp(0.15*log(1+(100/la)));
ra=0.856*exp(0.14*log(d));
re=1/(1+0.24*ec);

for x=1:n,
    for y=1:n,
        % calcul de la fréquence (cycles/pixel)
        x1=x-((n+1)/2);
        y1=y-((n+1)/2);
        f=sqrt(((x1)*(x1))+((y1)*(y1)))/((n-1)*sqrt(2)); 
         % calcul de la fréquence (cycles/degre)
        f=(n*f)/w;
        % calcul de l  orientation en radian
        if (x1==0)
            if y1>0 
                angle=pi/2;
            else
                angle=-pi/2;
            end
        else
        angle=atan(y1/x1);
        end
       rteta= 0.11*cos(4*angle)+0.89;     
        
   % calcul de S1(f,la,w)     
  f2=f;
  cs=f*f*la*la;
  if cs==0  cs=0;  else cs=exp(-0.3*log(cs)); end
  cs=3.23*cs;
  if cs==0  cs=0;  else cs=exp(5*log(cs)); end
  cs=cs+1;
  if cs==0  cs=0;  else cs=exp(log(cs)/5); end
  s1=cs*f2*al*0.9*exp(-0.9*bl*f2)*sqrt(1+0.06*exp(-0.9*bl*f2));
  
   % calcul de S1(f/r,la,w)  
  f2=f/(ra*re*rteta);
  cs=f2*f2*la*la;
  if cs==0  cs=0;  else cs=exp(-0.3*log(cs)); end
  cs=3.23*cs;
  if cs==0  cs=0;  else cs=exp(5*log(cs)); end
  cs=cs+1;
  if cs==0  cs=0;  else cs=exp(log(cs)/5); end
  s2=cs*f2*al*0.9*exp(-0.9*bl*f2)*sqrt(1+0.06*exp(-0.9*bl*f2));
  
  csf(x,y)=p*s1;
  if s2<s1 csf(x,y)=p*s2; end
             
    end
end

function [image1]=FiltrageCortex(imageorig,cortex,base,n)
% function [image]=FiltrageCortex(imageorig,cortex,base,n)
% imageorig: image ou matrice origine
% cortex: les 30 filtres cortex
% base: filtre bande de base.
% n: dimension de l image
% image: 31 images 1:base  2:(1,1)  3:(1,2)   n:(k,l) avec n: k=(n-1/6)+1 l=reste(n-1/6);
% Notes:les  valeurs de "image" sont réelles format double.
%
%
%
% voir aussi:
%           DecompositionCortex    Dalycortex    Dalycsf    DalyLuminance   Daly   ListeDaly   

csf=Dalycsf(n);
filtre=ones(n);
F=fftshift(fft2(imageorig));
F=F.*csf;
% filtrage bande de base
F3=F.*base;
F2=abs(ifft2(F3));
image1(1,:,:)=F2;

%for i=1:5,
%    for j=1:6,
%for x=1:n,
%    for y=1:n,
%  filtre(x,y)=cortex(i,j,x,y);
%    end
%end
%F3=F.*filtre;
%F2=abs(ifft2(F3));
%image1((i-1)*6+j+1,:,:)=F2;
%end
%end


for i=1:5
  for j=1:6
    filtre=reshape(cortex(i,j,:,:),256,256);
    F3=F.*filtre;
    F2=abs(ifft2(F3));
    image1((i-1)*6+j+1,:,:)=F2;
  end;
end;
