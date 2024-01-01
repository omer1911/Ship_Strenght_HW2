import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
offset = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.09, 0.324],
        [0.011, 0.055, 0.071, 0.083, 0.114, 0.335, 0.604],
        [0.047, 0.131, 0.173, 0.218, 0.314, 0.581, 0.813],
        [0.203, 0.343, 0.458, 0.573, 0.772, 0.98, 1.094],
        [0.439, 0.632, 0.809, 0.978, 1.089, 1.151, 1.187],
        [0.693, 0.908, 1.085, 1.18, 1.2, 1.2, 1.199],
        [0.856, 1.049, 1.183, 1.2, 1.2, 1.2, 1.2],
        [0.768, 0.993, 1.151, 1.188, 1.193, 1.191, 1.193],
        [0.502, 0.788, 0.986, 1.066, 1.087, 1.101, 1.123],
        [0.211, 0.484, 0.672, 0.751, 0.787, 0.82, 0.863],
        [0.041, 0.194, 0.295, 0.341, 0.357, 0.374, 0.429],
        [0.013, 0.084, 0.135, 0.149, 0.155, 0.165, 0.207],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.024]])*10     # hidro yarı genişlik değerleri 
boy = 120
genislik = 24
draft = 8.73 
yogunluk = 1.025
derinlik = 1.5 * draft
D=13.09
Cb=0.61
rho=1.025
posta0 = np.array([0, .5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5, 10]) * boy / 10
posta = np.linspace(0, boy, 201)
suhatti = np.array([0, .3, 1, 2, 3, 4, 5]) * draft / 4
w = boy * genislik * draft * Cb * rho # deplasman hesabı

def bonjeanAlani(offset, posta0, posta, suhatti):
    offset_new = np.zeros((201, 7))
    for i in range(7):
        f = interp1d(posta0, offset[:, i], kind = "cubic")
        offset_new[:, i] = f(posta)
    alan = np.zeros((201, 7))   # BON-JEAN ALANLARI
    for i in range(201):
        alan[i, 1:] = 2 * cumtrapz(offset_new[i, :], suhatti[:])
    return alan

def ataletDagilimi(boy, posta):
    Ix = np.zeros(201)   # ATALET MOMENT DAĞILIMI
    Wmin =4.032 #minimum midship section modulus [GL LOYD SECTİON 5C 2.1]
    Iy = 3 * Wmin * boy / 100  # ORTA KESİT ATALET MOMENTİ [m4]
    for i in range(201):
        if posta[i] <= boy / 20:
            Ix[i] = 5 * Iy * posta[i] / boy
        elif boy / 20 < posta[i] <= 7 * boy / 20:
            Ix[i] = .25 * Iy + (15 * Iy) * (posta[i] - boy / 20) / (6 * boy)
        elif 7 * boy / 20 < posta[i] <= 15 * boy / 20:
            Ix[i] = Iy
        elif 15 * boy / 20 < posta[i] <= 19 * boy / 20:
            Ix[i] = Iy - 2.5 * (Iy / boy) * (posta[i] - 15 * boy / 20)
        elif posta[i] > 19 * boy / 20:
            Ix[i] = .5 * Iy - 10 * (Iy / boy) * (posta[i] - 19 * boy / 20)
    return Ix

def celik_tekne(boy,genislik,derinlik,cb):
    N = boy * genislik * derinlik
    cs = (.21 - .026 * np.log10(N)) * (1 + .025 * (boy / derinlik - 12))
    G = cs * N * (1 + (2 / 3) * (cb - .7))
    qxs=prohaskaDagilimi(boy,G)
    return qxs,G
def alan_eq(suhatti,alan):
    polinomlar = np.zeros(201, dtype=object)  # Her bir postadaki alanlar için denklem oluşturulması
    degree = 6  # 6 tane su hatti ile yapıldığı için derece 6 olarak belrilendi
    for i in range(201):
        coefficients = np.polyfit(suhatti, alan[i], degree)
        polinomlar[i] = np.poly1d(coefficients)
    return polinomlar

def prohaskaDagilimi(boy, deplasman):
    qx = np.zeros(201)   # TOPLAM GEMİ AĞIRLIK DAĞILIMI (PROHASKA YÖNTEMİ)
    a = (.68 * deplasman / boy)
    b = 1.185 * deplasman / boy
    c = (.58 * deplasman / boy)
    n = 40
    posta=np.linspace(0,120,201)
    for i in range(201):
        if i < 68:
            qx[i] = a + posta[i] * (b - a) / 40   #uzunluğa bölündü oran-orantı yapıldı
        elif i > 135:
            qx[i] = b - n * (b - c) / 40  #uzunluğa bölündü oran-orantı yapıldı
            n -= 0.6
    qx[67:136]=b   
    qx[136:]=qx[136:][::-1]  #listeye n ilk 40 değeri ile girdiği için c 136. eleman olarak geldi bu yüzden liste ters çevrildi.
    return qx

def hesap_su(name,ax,qx,posta):
    px = ax-qx
    if name == "accident" :
            dpx0 = np.array([0, *cumtrapz(px,posta)])
            dpx = dpx0 - dpx0[-1] * posta / boy# LİNEER DÜZELTME (son değer en büyük değerin 0.03 unden küçük)
            ddpx0 = np.array([0, *cumtrapz(dpx,posta)]) # LİNEER DÜZELTME (son değer en büyük değerin 0.06 sından küçük)           
            ddpx = ddpx0 - ddpx0[-1] * posta / boy 
            Qx = dpx*boy/200
            Mx = ddpx*(boy/200)**2
    elif name =='launching':
        dpx = np.array([0, *cumtrapz(px,  posta)])
        dppx = np.array([0, *cumtrapz(dpx,  posta)])
        Qx = dpx*boy/200
        Mx = dppx*(boy/200)**2
    return Qx,Mx

def grafik(name,  Qx, Mx,Stress,vonmises):
        plt.figure(figsize=(10, 4))
        plt.title(name, fontweight="bold")
        plt.xlabel("Gemi Boyu [m]", fontweight="bold")
        plt.plot(  posta, Qx, posta, Mx,posta,Stress,posta,vonmises)
        plt.grid(color="red", linestyle="-", linewidth=.3)
        plt.legend([ "Q [ton]", "M/10 [ton.m]","Stress[Mpa]","vonmisses[Mpa]"], loc="best")
        plt.show() 
Ix=ataletDagilimi(boy,posta)
Ix -= .15 * Ix
alan=bonjeanAlani(offset,posta0,posta,suhatti)
polinomlar=alan_eq(suhatti,alan)
ax = alan[:, 5] * rho
deplasman_new = np.trapz(ax, posta)
qxs,G=celik_tekne(boy,genislik,derinlik,Cb)
alan[:, 5][85:130] = np.zeros(45)
ax = alan[:, 5] * rho
qx=prohaskaDagilimi(boy,deplasman_new)
LCG = np.sum(qx * posta) / np.sum(qx)
LCB = np.sum(alan[:, 5] * posta) / np.sum(alan[:, 5])
Qx,Mx=hesap_su("accident",ax,qx,posta)
LCG = np.sum(qxs * posta) / np.sum(qxs)
LCB = np.sum(alan[:, 5] * posta) / np.sum(alan[:, 5])
qx_n = np.trapz(qx, posta)
ymax = .6 * derinlik

A1 = derinlik * genislik
# Tarafsız eksenin 0.4xD old. kabul edildi
# Baş papet geminin tam başında 0.posta old. kabul edildi
# Kesme kuvvetin max old. postada max kayma gerilmesi
# Bu postada eninde perde olduğunu kabul edildi
# Perde sanki dikörtgen kesitli gibi hesaplama yapıldı
n = -10 # baş papatten bir kaç potsa önce
A2 = .4 * derinlik * genislik
S = (A1 - A2) * (.6 * derinlik) / 2
kayma = -9.81 * Qx[n] * S / (genislik * Ix[n] * 1000)
gerilme = 9.81 * Mx[1 : -1] * ymax / (Ix[1 : -1] * 1000)
gerilme = np.insert(gerilme, 0, 0)  
gerilme = np.append(gerilme, 0) 
von_mises = np.sqrt((gerilme**2) + (3 * kayma**2))
grafik("accident",Qx,Mx/10,gerilme,von_mises)
print("Maksimum Eğilme Gerilmesi =",max(abs(gerilme)))

M1=(0.9*boy-LCG)*G
print(M1)
egim= 0.044
ksi=np.zeros(201)
Batan_alan=np.zeros(201)
for i in range(169):
    ksi[i]=egim*(168*0.6-posta[i])
    Batan_alan[i]=polinomlar[i](ksi[i])

deplasman_q=Batan_alan*rho  #alanların oluşturduğu deplasman kuvvetleri ton/m
deplasman_new2=np.trapz(deplasman_q,posta)
LCB=np.sum(Batan_alan * posta) / np.sum(Batan_alan)
M2=(0.9*boy-LCB)*deplasman_new2
Qx,Mx=hesap_su("launching",deplasman_q,qxs,posta)
print(abs(max(Mx)))
gerilme = 9.81 * Mx[1 : -1] * ymax / (Ix[1 : -1] * 1000)
gerilme = np.insert(gerilme, 0, 0)  
gerilme = np.append(gerilme, 0) 
von_mises = np.sqrt((gerilme**2) + (3 * kayma**2))
print("Maksimum Eğilme Gerilmesi =",max(abs(gerilme)))
grafik("launching",Qx,Mx/10,gerilme,von_mises)







