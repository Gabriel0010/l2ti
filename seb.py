def debut():
  """import seb; plt,np,sig = seb.debut()"""  
  import matplotlib.pyplot as plt
  params = {'legend.fontsize': 20,
         'axes.labelsize': 20,
         'axes.titlesize':20,
         'xtick.labelsize':20,
         'ytick.labelsize':20,
         'legend.loc':'upper right'}
         
  plt.rcParams.update(params)
  #plt.rcParams['text.usetex'] = True
  import numpy as np
  #import scipy.signal as sig
  return plt,np

#################################################################################
#save and load
def save(nom_fichier,list_nom_var,list_var):
  """sauvegarde sous format binaire la liste des variables indiquees dans list_var
  le fichier s'appelle nom_fichier
  les noms des variables doivent etre mises avec des apostrophes autour
  """
  import pickle
  assert len(nom_fichier)>4, 'nom_fichier doit faire plus que 4 lettres'
  assert nom_fichier[-4:] == '.pkl'
  assert type(list_var) == list
  assert type(list_nom_var[0]) == str
  assert len(list_var) == len(list_nom_var)
  list_var_=[list_nom_var,list_var]
  open_file = open(nom_fichier, "wb")
  pickle.dump(list_var_, open_file)
  open_file.close()
  
def load(nom_fichier):
  """sauvegarde sous format binaire la liste des variables indiquees dans list_var
  le fichier s'appelle nom_fichier
  """
  import pickle
  assert len(nom_fichier)>4, 'nom_fichier doit faire plus que 4 lettres'
  assert nom_fichier[-4:] == '.pkl'
  file = open(nom_fichier, "rb")
  list_var = pickle.load(file)
  dc={}
  for k in range(len(list_var[0])):
    dc[list_var[0][k]]=list_var[1][k]
  file.close()
  return dc

def isfile(nom_fichier):
  """rend disponible os.path.isfile 
  True si nom_fichier existe
  """
  import os
  return os.path.isfile(nom_fichier)

########################################################################################
#Fonctions de bases
  
def fonction_T(t):
  """y=fonction_T(t)"""
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray or np.isscalar(t)
  y=(1-np.abs(t))*(np.abs(t)<=1)
  return y


def fonction_H(t):
  """y=fonction_H(t) echelon"""
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray or np.isscalar(t)
  y=(t>=0)
  return y

  
def fonction_P(t):
  """y=fonction_P(t) porte
    =[|t|<=0.5](t)   
  """
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray or np.isscalar(t)
  y=(np.abs(t)<=0.5)
  return y
  
  
def fonction_C(t):
  """y=fonction_C(t) 
  """
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray or np.isscalar(t)
  y=(0.5+t)*((t>=-1/2)&(t<=1/2))
  return y
  
def fonction_D(t):
  """y=fonction_D(t) 
  """
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray or np.isscalar(t)
  y=(0.5-t)*((t>=-1/2)&(t<=1/2))
  return y

######################################################################################
#Filtrage
def convolution(tx,x,th,h,ty):
  """Le programme fournit le produit de convolution de 
     x par h aux instants demandes par ty
     tx, th et ty sont les echelles de temps de x h et y
     x et tx doivent etre de meme taille
     th et h doivent etre de meme taille
     tx doit contenir au moins deux composantes
     x et h sont supposés d'énergie finie. 
  """
  import seb
  plt,np = seb.debut()
  #import scipy.signal as sig
  assert len(tx)==len(x), 'tx doit avoir une echelle de temps compatible avec x '
  assert len(th)==len(h)
  Te=tx[1]-tx[0]   
  assert _regulierement_reparti_(tx)
  if len(ty)>1:
    assert _est_Te_(ty,Te)
  if len(th)>1: 
    assert _est_Te_(th,Te)
  dtx = _delta_t_(tx,Te); tx = tx-dtx
  dty = _delta_t_(ty,Te); ty = ty-dty
  dth = _delta_t_(th,Te); th = th-dth
  assert np.abs(_delta_t_(tx,Te))<1e-10
  assert np.abs(_delta_t_(ty,Te))<1e-10 
  assert np.abs(_delta_t_(th,Te))<1e-10
  yconv = np.convolve(x,h)
  assert len(yconv) == len(x) + len(h) - 1
  y = np.zeros_like(ty)
  tconv = np.arange(tx[0]+th[0],tx[-1]+th[-1]+1e-10,Te)
  assert len(yconv) == len(tconv) 
  #signal y avant yconv
  if max(ty)<tconv[0]-1e-10:
    return y
  #signal y au milieu de yconv
  elif (min(ty)>tconv[0]+1e-10) and (max(ty)<=tconv[-1]+1e-10):
    idebut = np.nonzero(tconv >= min(ty)-1e-10)[0][0]
    iapres_fin = len(ty)+idebut 
    assert iapres_fin-idebut == len(ty)
    y=yconv[idebut:iapres_fin]
  #signal y a cheval sur la fin    
  elif (min(ty)>tconv[0]+1e-10) and (max(ty)>tconv[-1]+1e-10):
    idebut = np.nonzero(tconv >= min(ty)-1e-10)[0][0]
    iapres_fin=len(tconv)+2-idebut
    y[:iapres_fin]=yconv[idebut:]
  #signal y a cheval sur le debut
  elif max(ty)<=max(tconv)+1e-10:
    idebut=np.nonzero(ty>=min(tconv)-1e-10)[0][0]
    iapres_fin=len(ty)-idebut
    y[idebut:]=yconv[:iapres_fin]
  #signal y autour   
  else: 
    assert min(ty)<min(tconv)
    assert max(ty)>max(tconv)
    idebut=np.nonzero(ty>=min(tconv)-1e-10)[0][0]
    assert ty[idebut]>=min(tconv)-1e-10
    iapres_fin=len(tconv)+idebut
    assert ty[iapres_fin-1]>=max(tconv)-1e-10
    y[idebut:iapres_fin]=yconv
  y=y*Te
  return y

def correlation(tx,x,ty,y,tz):
  """calcule l'intercorrelation entre tx,x et ty,y en tz"""
  import seb
  plt,np = seb.debut()
  assert len(ty.shape)==1 ,('le vecteur ty doit etre colonne c-a-d une dimension',ty.shape,len(ty.shape))
  assert len(y.shape)==1 ,('le vecteur ty doit etre colonne c-a-d une dimension',y.shape,len(y.shape))
  ty2=np.flipud(-ty)
  y2=np.flipud(y)
  gamma=seb.convolution(tx,x,ty2,y2,tz)
  return gamma


#########################################################################################
#Equation differentielle 

def sol_eq_diff(coef,t):
  """y=sol_eq_diff(coef,t)
  coef sont les coefficients devant les termes de l'equation differentielle definis 
  comme un tuple
  t est l'ensemble des instants dont a cherche a calculer y(t)
  t est un vecteur avec des points regulierement espaces
  """
  assert type(coef) == tuple
  import seb
  plt,np = seb.debut()
  assert type(t) == np.ndarray
  t1,t2 = _cut_t_(t)
  assert len(t1)+len(t2) == len(t)
  assert max(t1) <= 0
  assert min(t2) >= 0
  assert _is_in_(t1,t)
  assert _is_in_(t2,t)
  import scipy.signal as sig
  y1 = np.zeros_like(t1)
  if 0 == len(t2):
    y2 = np.array([])
  else:
    t3,y2,_ = sig.lsim(((1,0),coef),np.ones_like(t2),t2)
  assert len(t3) == len(t2)  
  y=np.concatenate((y1,y2))
  return y
  

def deriver(t,x):
  """derive le signal x(t)"""
  Te=t[1]-t[0]
  import numpy as np
  y=np.concatenate((np.diff(x)/Te,np.array([0])))
  return y

def integrer(t,x):
  """derive le signal x(t)"""
  Te=t[1]-t[0]
  import numpy as np
  y=np.cumsum(x)*Te
  return y

############################################################################
#Transformer les signaux

def retarder(t,x,tau):
  """Cette fonction retarde le signal x(t) defini par t et x
  de tau lorsque tau est positif et avance de -tau si tau est negatif
  """
  Te=t[1]-t[0]
  import numpy as np
  assert np.isscalar(tau), tau
  if tau<0:
    ind=round(-tau/Te)
    y1=x[ind:]
    y2=np.zeros(ind)
    y=np.concatenate((y1,y2))
    # print(len(y),len(x),ind,len(y1),len(y2))
    assert len(y)==len(x)
  elif tau>0:
    ind=round(tau/Te)
    # print(type(ind),ind)
    y1=x[:len(x)-ind]
    y2=np.zeros(ind)
    y=np.concatenate((y2,y1))
    # print(len(y),len(x),ind,len(y1),len(y2))
    assert len(y)==len(x)
  else:
    y = x
  return y    


def periodiser_ech_t(t,T):
  """t2=periodiser_ech_t(t1,T)
    produit un vecteur de meme taille que t1 mais dont les 
    valeurs sont entre 0 et T de facon a definir un 
    signal periodique de periode T
    Si T est un intervalle alors c'est le motif entre T[0] et T[1] qui est repete
  """
  import seb
  (plt,np) = seb.debut()
  assert type(t) == np.ndarray, (t,type(t))
  if np.isscalar(T):
    return t-np.floor(t/T)*T
  else: 
    assert type(T) == tuple, ('Si T n est pas la periode alors c est un tuple decrivant l intervalle ',T)
    assert len(T)  == 2,     ('Si T n est pas la periode alors c est un doit avoir deux composantes',T)
    P=T[1]-T[0]
    return t-np.floor((t-T[0])/P)*P

#def echantillonner(tx,x,fe):
#  """echantillonne le signal defini par t,x a la frequence fe
#  le signal echantillonne est defini par te,xe
#  """
#  Te=1/fe

################################################
#outils generaux 
def find_nearest(array, value):
  import numpy as np
  array = np.asarray(array)
  array = array[np.isfinite(array)]
  idx = (np.abs(array - value)).argmin()
  return array[idx]

def where_nearest(array, value):
  import numpy as np
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return idx

def gaussian(x,mu,sigma):
  """renvoie la fonction associe a densite de probabilite de la gaussienne de moyenne mu et d ecart type sigma
  """
  from scipy.stats import norm
  rv = norm(loc = mu, scale = sigma)
  return rv.pdf(x)

def val(x):
  """vérifie si x est un numpy array contenant une seule valeur, 
  une seule valeur ou autre chose
  Si c'est autre chose, cela met une erreur.
  Sinon cela renvoie cette unique valeur
  """
  import numpy as np
  if np.isscalar(x):
    return x
  elif type(x)==np.array and 1==len(x):
    return x[0]
  else:
    try:
      print(f"x est de type {type(x)}, de longueur {len(x)}")
    except:
      assert False
      
####################################################################
#Calcul erreur quadratique   
def erreur_quad(fun,intervalle):
  """cette fonction utilise une variable aleatoire sur un support uniforme pour 
  calculer une erreur quadratique. 
  intervalle doit indique avec un tuple contenant deux valeurs
  """
  import seb
  plt,np = seb.debut()
  assert type(intervalle) == tuple
  assert len(intervalle) == 2
  T=np.random.uniform(intervalle[0],intervalle[1],10**5)
  assert np.isreal(fun(T[0]))
  assert fun(T[0])**2 >= 0
  return np.sqrt(np.sum(fun(T)**2))

  
########################################################################### 
#Transformees de Fourier 

def TFTD(t,x,f):
  """TFTD du signal temps discret t,x en les fréquences f
  """
  import numpy as np
  Te=t[1]-t[0]
  assert np.abs(t[1]/Te-np.round(t[1]/Te))<1e-8, ("t doit etre un multiple de la periode d'echantillonnage",'t[1]/(t[1]-t[0])=',t[1]/(t[1]-t[0]))
  if np.isscalar(f):
    X=np.sum(x*np.exp(-1j*2*np.pi*f*t))
  else:
    X=np.zeros(len(f)).astype(complex)
    for f_ in range(len(f)):
      X[f_]=np.sum(x*np.exp(-1j*2*np.pi*f[f_]*t))

  return X  

def TFTDI(f,X,t):
  assert len(X.shape) == 1, "x ne doit pas etre une matrice"
  assert f.shape == X.shape, ("f et X doivent avoir la meme taille",f.shape,X.shape)
  import seb
  plt,np = seb.debut()
  if 0==len(t):
    return np.array([])
  elif np.isscalar(t):
    fe=f[-1]-f[0]
    x=np.trapezoid(X*np.exp(1j*2*np.pi*f*t),f)/fe
  else:
    assert t.shape==(len(t),) ,(t.shape)
    assert len(t)>0, len(t)
    fe=1/(t[1]-t[0])
    assert not np.isscalar(f)
    assert np.abs(f[-1]-f[0] - fe)<fe/len(f)*10
    x=np.zeros(t.shape,dtype=np.complex128)
    for t_ in range(len(t)):
      x[t_] = np.trapezoid(X*np.exp(1j*2*np.pi*f*t[t_]),f)/fe
  return x

  
def TFI(f,X,t):
  """x=TF(f,X,t) est la transformee de Fourier inverse  du spectre defini par f et X en t
  """  
  assert len(X.shape) == 1, ("X ne doit pas etre une matrice",X.shape,len(X.shape))
  assert f.shape == X.shape, "t et x doivent avoir la meme taille"
  import seb
  x=seb.TF(f,X,-t)
  return x  



def TF(t,x,f):
  """X=TF(t,x) est la transformee de Fourier du signal defini par t et x
  """
  assert t[1] >= t[0], "t doit etre dans un ordre croissant"
  assert len(x.shape) == 1, "x ne doit pas etre une matrice"+str(x)
  assert t.shape == x.shape, "t et x doivent avoir la meme taille"
  import seb
  plt,np = seb.debut()
  X=seb.TF1(t,x,f)
  assert t[1]-t[0]>=10**(-8),   ('t=',t)
  fe=1/(t[1]-t[0])
  assert np.abs(np.sinc(1))< 1e-13
  X=X*np.sinc(f/fe)*np.exp(-1j*np.pi*f/fe) #X=X*np.sinc(f/fe)
  return X

def TF1(t,x,f):
  """TF classique TFTD et methode de trapeze et multiplilcation par Te
  utilise par seb.TF
  """
  assert len(x.shape) == 1, "x ne doit pas etre une matrice"
  assert t.shape == x.shape, "t et x doivent avoir la meme taille"
  import seb
  plt,np = seb.debut()
  if np.isscalar(f):
    X=np.trapezoid(x*np.exp(-1j*2*np.pi*f*t),t)
  else:
    X=np.zeros(f.shape,dtype=np.complex128)
    for f_ in range(len(f)):
      X[f_] = np.trapezoid(x*np.exp(-1j*2*np.pi*f[f_]*t),t)
  return X
  
def argmax_fdf(f,X):
  """calcule la frequence maximisant le spectre et le plus petit décalage permettant de diviser par deux le module 
  (f,delta_f)=seb.argmax_fdf(f,X)
  """
  import seb
  plt,np = seb.debut()  
  f0_=np.argmax(abs(X))
  f0=abs(f[f0_])
  Xp=X[f0_:-1]
  Xm=X[::-1][f0_:-1]
  l=min(len(Xp),len(Xm))
  Xp=Xp[1:l]
  Xm=Xm[1:l]
  assert Xp.shape == Xm.shape, "Xp et Xm doivent avoir la meme taille"
  test=np.logical_or(abs(Xm)<=0.5*abs(X[f0_]),abs(Xp)<=0.5*abs(X[f0_]))
  if any(test):
    delta_f_=np.argwhere(test)[0][0]
  else: 
    delta_f_=l
  delta_f=f[f0_+delta_f_]-f[f0_]
  assert delta_f>=0
  return f0,delta_f

def coef_serie_Fourier(t,x,T,k):
  """calcul les coefficients de la serie de Fourier Xk
  T est soit la periode du signal soit un tuple indiquant un intervalle sur lequel est defini x(t)
  Si T est une valeur alors l'intervalle considere est [0,T]
  k est la liste des indices des frequences calcules
  Le programme retourne un tuple avec d'abord les frequences et d'autre part les coefficients
  """
  # assert False, "cas chevauchement non traité"
  import numpy as np
  import seb
  if np.isscalar(T): 
    assert T>0,                           ('T doit etre positif', T)    
    T=(0,T); P=T[1]-T[0]
  else: 
    assert 2==len(T)    
    P=T[1]-T[0]
  assert all(np.abs(np.round(k)-k))<1e-10, ('k doit etre une liste d entiers',k)
  t1=t[(t>=T[0])&(t<=T[1])]
  x1=x[(t>=T[0])&(t<=T[1])]
  f=k/P
  if np.isscalar(f):
    X=np.trapezoid(x1*np.exp(-1j*2*np.pi*f*t1),t1)/T
  else:
    X=np.zeros(f.shape,dtype=np.complex128)
    for f_ in range(len(f)):
      X[f_] = np.trapezoid(x1*np.exp(-1j*2*np.pi*f[f_]*t1),t1)/P
  return f,X

    

def TFD(t,s,T,bool):
  """
  calcule la TFD, T indique soit la periode soit l'intervalle utilise pour decrire le signal periodique
  t,s defini l echelle de temps 
  bool vaut True si on veut une representation centree et False si on n'en veut pas 
  f est l echelle de frequence retournee 
  S est l ensemble des coefficients associes
  f,S=seb.TFD(t,s,T,bool)
  """
  import numpy as np
  import seb
  if np.isscalar(T): 
    assert T>0,                           ('T doit etre positif', T)    
    T=(0,T); P=T[1]-T[0]
  else: 
    assert 2==len(T)    
    P=T[1]-T[0]
  t1=t[(t>=T[0])&(t<=T[1])]
  s1=s[(t>=T[0])&(t<=T[1])]
  fe=1/(t[1]-t[0])
  assert np.abs(P*fe - np.round(P*fe))<1e-8, ('N doit etre entier',P,fe,np.abs(P*fe - np.round(P*fe)))
  N=round(P*fe)
  if 1==N:
    f=np.array([0]); 
    S=np.array([s1[0]])
    return f,S
  fe=N/P
  if bool:
    if N/2==round(N/2):
      f=np.concatenate((np.arange(-fe/2,0,fe/N),np.arange(0,fe/2,fe/N)))
    else: 
      f=np.concatenate((np.arange(-fe/2+fe/N/2,0,fe/N),np.arange(0,fe/2,fe/N)))
  else:
    np.arange(0,fe,fe/N)
  assert len(f)==N,         (len(f),N)
  assert np.abs(f[0]-seb.synchroniser(f)[0])<1e-8, (f,np.abs(f[0]-seb.synchroniser(f)[0]))
  assert max(f)<fe
  assert min(f)>=-fe/2
  assert all(np.abs(np.diff(f)-fe/N)<1e-8), (f,fe,N)
  import scipy.fft
  S_fft=scipy.fft.fft(s)/N
  if bool: 
    S=np.concatenate((S_fft[round(N/2):],S_fft[0:round(N/2)]))
  else: 
    S=S_fft
  assert len(f)==len(S)  
  return f,S

#####################################################################
#sur l'echelle de temps 
def synchroniser(t):
  """
  change l'echelle de temps de façon que le vecteur soit un multiple 
  de la periode d'echantillonnage
  """
  import numpy as np
  assert len(t)>=2, ('t doit avoir au moins deux composantes',len(t))
  Te=t[1]-t[0]; 
  return t+np.round(t[0]/Te)*Te-t[0]
  
########################################################################################## 
#Sous-programmes
    
def _est_Te_(t,Te):
  """private
  verifie que Te est bien la periode d'echantillonnage de t
  """
  import numpy as np
  return all(np.abs(np.diff(t) - Te)<1e-12)
  
def _regulierement_reparti_(t):
  """private
  verifie que le vecteur t est regulierement reparti
  """  
  import numpy as np
  return all(np.abs(np.diff(t) - (t[1]-t[0]))<1e-12)
  
def _cut_t_(t):
  """private
  decoupe le vecteur t en 
  deux vecteurs, le premier ayant les valeurs negatives et le deuxieme les valeurs positives. 
  """
  import numpy as np
  assert _regulierement_reparti_(t)
  if (max(t) >= 0) and (min(t) < 0): 
    ind = np.where(t>=0)[0][0]
    t2 = t[ind:]
    t1 = t[:ind]
  elif min(t) >= 0:
    t1 = np.array([])
    t2 = t
  else: 
    t1 = t 
    t2 = np.array([])
  return t1,t2

def _is_in_(a,b):
  """private 
  verifie que a est inclu dans b
  """
  assert len(a)>0
  assert len(b)>0
  import numpy as np
  assert type(a) == np.ndarray
  import random
  if 0 == len(a):
    return true
  else:
    a_ = random.choice(a)
    return len(np.where(b==a_))>0


def _delta_t_(t,Te):
  """private
  calcule le décalage entre 0 et le premier instant après 0
  """
  import numpy as np
  return t[0]-Te*np.floor(t[0]/Te+1e-13)

  
if __name__ == '__main__':
    debut()  
