function [alpha,lambda,LambdaBis,coord,intensite]=qualityMetric2Alpha(im,structure,measureName)
%    This program illustrates the following paper 
%    G. Dauphin and P. Viaris de Lesegno, Assessment of Quality Metrics with 
%    Reference Based on Uniform Colour Spaces. In Proceeding of 2nd European 
%    Workshop on Visual Information Processing. 6 pages, Paris, France, 
%    July 2010. EUVIP'10.
%
%    With this program it is possible to extract a parameter alpha from any
%    quality metrics, if the extracted parameter is close to 0.33 the quality 
%    metric under study is consistent with HVS regarding its sensitivity 
%    to background luminance. If the extracted parameter is close to 0 the 
%    quality metric under study is consistent with Weber's law. 
%
%    measureName has to be the name of a Matlab function with two parameters,
%    the first parameter is for the original image,
%    the second parameter is for the distorted image,
%    it should accept images defined as matrices of double,
%    it should be in the search path (use addpath to add a directory).
%
%    The original image should be a grey valued image defined 
%    as a matrix of double.
%
%    im=double(imread('fileName'))/256;
%
%    The structure is any binary matrix of fairly small size.
%
%    structure=[1 0 0 0 0 0 0;...
%               1 1 0 0 0 0 0;...
%               1 1 1 0 0 0 0;...
%               1 1 1 1 0 0 0;...
%               1 1 1 1 1 0 0;...
%               1 1 1 1 1 1 0;...
%               1 1 1 1 1 1 1 ];
%
%    structure=[-ones(1,8) ones(1,8)];
%
%    structure=ones(7);
%
%    The output of this function is the parameter alpha extracted, 
%    the relationship between lambda and lambda prime, 
%    the coordinate where the distortion has been embedded, 
%    the intensite in cd.m-2 with which it has been embedded.  
%
%    [alpha,lambda,LambdaBis,coord,intensite]=qualityMetric2Alpha(im,structure,measureName)
%
if ~(isImage(im)) 
  disp('The first parameter is not correct.'), 
  alpha=NaN;
  return; 
end;

coord=choseCoord(size(im));
if ~(isCoord(im,coord,structure)) 
  disp('The chosen coordinate are not valid.'), 
  alpha=NaN;
  return; 
end;

distortion=structure2distortion(structure,coord,size(im));

intensite=ChoixIntensite(mean2(im(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)));

[QualRef,monotonie]=CalculRef(measureName,im,distortion,intensite);
if isnan(QualRef) alpha=NaN; return; end;

[LambdaBis,lambda,VectQual]=lambda2LambdaBis(measureName,im,distortion,intensite,QualRef,monotonie);

[alpha,beta]=CalculAlpha(lambda,LambdaBis);

function test=isImage(im)
  test=((2==length(size(im)))&(max(max(im))<3));

function test=isCoord(im,coord,structure)
  test=((coord(1)>=1+size(structure,1)/2)&(coord(1)<=size(im,1)-size(structure,1)/2)&(coord(2)>=1+size(structure,2)/2)&(coord(2)<=size(im,2)-size(structure,2)/2));

function coord=choseCoord(dimension)
%this function decides where the structure is being embedded.
  coord(1)=ceil(rand(1)*dimension(1));
  coord(2)=ceil(rand(1)*dimension(2));

function intensite=ChoixIntensite(gray)
  lum2gray=inline('sign(L).*(abs(L)/180).^(1/2.4)','L');
  intensite=lum2gray(rand(1)*3.5+0.5+gray)-lum2gray(gray);

function distortion=structure2distortion(structure,coord,dimension)
  distortion=zeros(dimension);
  i1=coord(1)-floor(size(structure,1)/2);
  i2=i1+size(structure,1)-1;
  j1=coord(2)-floor(size(structure,2)/2);
  j2=j1+size(structure,2)-1;
  distortion(i1:i2,j1:j2)=structure;

function [QualRef,monotonie]=CalculRef(NomMesure,im,distortion,intensite)
  gray2lum=inline('180*abs(g).^(2.4).*sign(g)','g'); % gamma=2.4;L0=180;
  lum2gray=inline('sign(L).*(abs(L)/180).^(1/2.4)','L');
  ImRef=im;
  ImDis=lum2gray(gray2lum(im)+intensite*distortion);
  try
    eval(['QualRef=',NomMesure,'(ImRef,ImDis);']);
  catch
    disp('The third parameter seems to be wrong or has not been found'),
    QualRef=NaN; monotonie=NaN;
  end;
  ImDis=lum2gray(gray2lum(im)+intensite/2*distortion);
  try
    eval(['QualBis=',NomMesure,'(ImRef,ImDis);']);
  catch
    disp('The third parameter seems to be wrong or has not been found'),
    QualRef=NaN; monotonie=NaN;
  end;
  monotonie=sign(QualRef-QualBis);

function [LambdaBis,lambda,VectQual]=lambda2LambdaBis(NomMesure,im,distortion,intensite,QualRef,monotonie)
  Dichotomy=inline('[LumSet(1) 0.5*LumSet(1)+0.5*LumSet(2)]*test+[0.5*LumSet(1)+0.5*LumSet(2) LumSet(2)]*(1-test)','LumSet','test');
  gray2lum=inline('180*abs(g).^(2.4).*sign(g)','g'); % gamma=2.4;L0=180;
  lum2gray=inline('sign(L).*(abs(L)/180).^(1/2.4)','L');
  lambda=(0.5:10.5)/11*gray2lum(1)/max(max(gray2lum(im)));
  LambdaInit=[0 5];
  for k=1:length(lambda)
    LambdaSet=LambdaInit;
    for l=1:15
      ImRef=lum2gray(lambda(k)*gray2lum(im));
      ImDis=lum2gray(lambda(k)*gray2lum(im)+...
	mean(LambdaSet)*intensite*distortion);
      try
        eval(['qual=',NomMesure,'(ImRef,ImDis);']); 
      catch
        disp('The third parameter seems to be wrong or has not been found'),
        QualRef=NaN; monotonie=NaN;
      end;
      if (1==monotonie) 
        LambdaSet=Dichotomy(LambdaSet,(qual>QualRef));
      else LambdaSet=Dichotomy(LambdaSet,(qual<QualRef));
      end;
    end;
    VectQual(k)=qual;
    LambdaBis(k)=mean(LambdaSet);
  end;

function [alpha,beta]=CalculAlpha(lambda,LambdaBis)
  LambdaInit=[0 5];
  SetShorter=inline('[u(1)*(1-0.5*r)+u(2)*0.5*r u(1)*0.5*r+u(2)*(1-0.5*r)]','u','r');
  FindInSet=inline('find((X>min(Set))&(X<max(Set)))','X','Set');
  CorAf=inline('[sum(x(:).*y(:))*length(x)-sum(y)*sum(x) sum(x.^2)*sum(y)-sum(x(:).*y(:))*sum(x)]/(length(x)*sum(x.^2)-sum(x)^2)','x','y');
  index=FindInSet(LambdaBis(:),SetShorter([min(LambdaInit) max(LambdaInit)],1/2^14));
  vect=CorAf(log(lambda(index)+eps),log(LambdaBis(index)+eps));
  alpha=1-vect(1); 
  beta=vect(2);

