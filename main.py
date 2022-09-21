import streamlit as st
import pandas as pd
# from functions import *
from module import *
import matplotlib.pyplot as plt
from PIL import Image
# plt.style.use('science')
import io
from pathlib import Path, PurePath
from numba import njit




image1 = Image.open('cismid.jpg')
st.image(image1)
st.write('# Escalamiento de Espectros')

image2 = Image.open('uni.jpg')
st.sidebar.image(image2,width=150)
st.sidebar.header('## DATOS')
st.sidebar.markdown("""
[Ejemplo de Archivo de Entrada](https://github.com/elvis1398/SpectrumScal/blob/main/PRQ_19661017164100.txt)
""")

if st.button('Borrar Registros Anteriores'):
    try:
        f = open('.\\Espectros.txt', "w")
        f.write('')
        a = open('.\\Acelerogramas.txt', "w")
        a.write('')
    except:
        pass
    

# Cargar Archivos
uploaded_files = st.sidebar.file_uploader("Choose a TXT file",accept_multiple_files=True)
# st.write(uploaded_files)

sismo = []
nombre_registros = []
for archivo in uploaded_files:
    data = np.loadtxt(archivo,skiprows=37)
    sismo = sismo + list(data)
    nombre_registros += [archivo.name]
    
sismo = sismo + [[0.0,0,0,0]]
sismo = np.array(sismo)
    
# np.savetxt('.\\Registros.txt',sismo)

    
nregistros = len(np.argwhere(sismo[:,0]==0.0))-1
i0 = np.zeros(nregistros+1)
# st.write(i0)

for i in range(nregistros+1):
    i0[i] = int(round(np.argwhere(sismo[:,0]==0.0)[i][0],0))
# st.write(i0)


# datos = pd.DataFrame(data, columns=['t','EO','NS','UD'])
# st.write('Datos',datos)

# Parámetros
with st.sidebar.container():
    Test = st.number_input('Periodo de la Estructura', 0.05, 1.5, 0.2)
    
    
direccion = st.sidebar.radio("Dirección del Sismo",('EO', 'NS'))

# if direccion == 'EO':
#     Ag0 = data[:,1]
# else:
#     Ag0 = data[:,2]
# try:
if uploaded_files[0] is not None:
    if open(".\\Acelerogramas.txt",'r').read()=='':    
        #Frecuencia de Nyquist, discretización
        Ag0_EO = np.zeros(sismo.shape[0])
        Ag0_NS = np.zeros(sismo.shape[0])
        dt0 = np.zeros(nregistros)
        ndatos0 = np.zeros(nregistros)

        Ag_EO = []
        Ag_NS = []
        dt = []
        ndatos = []

        t = []

        for i in range(nregistros):
            Ag0_EO[int(i0[i]):int(i0[i+1])] = sismo[int(i0[i]):int(i0[i+1]),1]
            Ag0_NS[int(i0[i]):int(i0[i+1])] = sismo[int(i0[i]):int(i0[i+1]),2]
            dt0[i] = sismo[int(i0[i]+1),0] - sismo[int(i0[i]),0]
            ndatos0[i] = len(Ag0_EO[int(i0[i]):int(i0[i+1])])

            aux_EO,aux_ndatos,aux_dt = intsignal(dt0[i], Ag0_EO[int(i0[i]):int(i0[i+1])], int(ndatos0[i]))
            aux_NS,aux_ndatos,aux_dt = intsignal(dt0[i], Ag0_NS[int(i0[i]):int(i0[i+1])], int(ndatos0[i]))

            aux_t = np.zeros(aux_ndatos)
            for i in range(aux_ndatos):
                aux_t[i] = i*aux_dt


            Ag_EO = Ag_EO + list(aux_EO)
            Ag_NS = Ag_NS + list(aux_NS)
            dt = dt + [aux_dt]
            ndatos = ndatos + [aux_ndatos]
            t = t + list(aux_t)
        t = t +[0]

        Ag_EO = np.array(Ag_EO)
        Ag_NS = np.array(Ag_NS)
        dt = np.array(dt)
        ndatos = np.array(ndatos)    

        t = np.array(t)

        i0 = np.zeros(nregistros+1)
        # st.write(np.argwhere(t==0))

        for i in range(nregistros+1):
            i0[i] = int(round(np.argwhere(t==0)[i][0],0))
        # st.write(i0)

        #Filtrado y correción por línea base
        Ag_corr_EO = []
        Ag_corr_NS = []

        for i in range(nregistros):
            #Frecuencia de Nyquist, discretización
            aux_EO = Filt_Corr(Ag_EO[int(i0[i]):int(i0[i+1])],dt[i],ndatos[i],f1=0.1,f2=25)
            aux_NS = Filt_Corr(Ag_NS[int(i0[i]):int(i0[i+1])],dt[i],ndatos[i],f1=0.1,f2=25)

            Ag_corr_EO = Ag_corr_EO + list(aux_EO)
            Ag_corr_NS = Ag_corr_NS + list(aux_NS)

        Ag_corr_EO = np.array(Ag_corr_EO)
        Ag_corr_NS = np.array(Ag_corr_NS)

        np.savetxt('.\\Acelerogramas.txt',
                   np.transpose((np.append(t[:-1],np.append(np.append(Ag_EO,Ag_NS),np.append(Ag_corr_EO,Ag_corr_NS)))).reshape((5,len(Ag_corr_NS)))))

    else:
        t, Ag_EO, Ag_NS, Ag_corr_EO, Ag_corr_NS = carchivos('.\\Acelerogramas.txt')
        # t = np.loadtxt('.\\Acelerogramas.txt')[:,0]
        # t = np.append(t,np.array([0.0]))
        # Ag_EO = np.loadtxt('.\\Acelerogramas.txt')[:,1]
        # Ag_NS = np.loadtxt('.\\Acelerogramas.txt')[:,2]
        # Ag_corr_EO = np.loadtxt('.\\Acelerogramas.txt')[:,3]
        # Ag_corr_NS = np.loadtxt('.\\Acelerogramas.txt')[:,4]

        ndatos = np.zeros(nregistros)
        dt = np.zeros(nregistros)

        i0 = np.zeros(nregistros+1)
        # st.write(np.argwhere(t==0))

        for i in range(nregistros+1):
            i0[i] = int(round(np.argwhere(t==0.0)[i][0],0))

        for i in range(nregistros):
            ndatos[i] = len(Ag_EO[int(i0[i]):int(i0[i+1])])
            dt[i] = t[int(i0[i]+1)]-t[int(i0[i])]




if direccion == 'EO':
    Ag = Ag_EO
    Ag_corr = Ag_corr_EO
else:
    Ag = Ag_NS
    Ag_corr = Ag_corr_NS

#Espectro de la norma de diseño Sismorresistente E030
periodo_P = {'S0':0.3,'S1':0.4,'S2':0.6,'S3':1.0,} #TP
periodo_L = {'S0':3.0,'S1':2.5,'S2':2.0,'S3':1.6,} #TL
suelo  = {'S0':{'Z4':0.8,'Z3':0.8,'Z2':0.8,'Z1':0.8,},
          'S1':{'Z4':1.0,'Z3':1.0,'Z2':1.0,'Z1':1.0,},
          'S2':{'Z4':1.05,'Z3':1.15,'Z2':1.20,'Z1':1.60,},
          'S3':{'Z4':1.10,'Z3':1.20,'Z2':1.40,'Z1':2.00,},}
zona = {'Z1':0.10,'Z2':0.25,'Z3':0.35,'Z4':0.45,}


with st.sidebar.container():
    st.write("Parámetros Sísmicos")
    col1, col2, col3 = st.columns(3)
    with col1:
        Z = st.selectbox('Z',('Z1', 'Z2','Z3', 'Z4'),index=3)
    with col2:
        u = st.selectbox('U',(1, 1.3, 1.5),index=0)
    with col3:    
        S = st.selectbox('S',('S0','S1','S2','S3'),index=0)

s = suelo[S][Z]
z = zona[Z]
TP = periodo_P[S]
TL = periodo_L[S]

T = np.arange(0.0,5.01,0.01,dtype = 'double')
T[0] = 0.02
C = np.zeros(T.shape[0])

ii = 0

for i1 in T:

    if i1 < 0.2*TP:
        C[ii] = 1+7.5*i1/TP
    elif 0.2*TP <= i1 <TP:
        C[ii] = 2.5
    elif TP <= i1 <TL:
        C[ii] = 2.5*TP/i1
    else:
        C[ii] = 2.5*TP*TL/i1**2
    ii += 1 

Sa_e030 = z*u*s*C*981 #gal



#Gráficos
fig1, ax1 = plt.subplots(nregistros,1,constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(15,4*nregistros))
fig2, ax2 = plt.subplots(nregistros,1,constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(15,4*nregistros))
fig3, ax3 = plt.subplots(nregistros,1,constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(15,8*nregistros))
fig4, ax4 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(15,8))
# fig5, ax5 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))    

tab1, tab2, tab3, tab4 = st.tabs(["Señales Originales", "Señales Filtradas y Corregidas","Espectros: SRSS", "Espectro Promedio"])



with tab1:
    for i in range(nregistros):
        ax1[i].plot(t[int(i0[i]):int(i0[i+1])],Ag[int(i0[i]):int(i0[i+1])], 'black', linewidth=1) #,label='original')
        # ax1.plot(t,Ag_corr, 'blue', linewidth=1,label='corregida')
        ax1[i].set_title(f'Aceleración del Suelo {nombre_registros[i]}', fontsize=20, fontweight = 'bold')
        ax1[i].set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
        ax1[i].set_ylabel('Aceleración (gal)', fontsize=15, fontweight = 'bold')
        ax1[i].tick_params(labelsize=12)
        ax1[i].legend(loc='best', fontsize=15)
    st.pyplot(fig1)



with tab2:
    for i in range(nregistros):
        ax2[i].plot(t[int(i0[i]):int(i0[i+1])],Ag_corr[int(i0[i]):int(i0[i+1])], 'blue', linewidth=1)
        ax2[i].set_title(f'Aceleración del Suelo {nombre_registros[i]}', fontsize=20, fontweight = 'bold')
        ax2[i].set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
        ax2[i].set_ylabel('Aceleración (gal)', fontsize=15, fontweight = 'bold')
        ax2[i].tick_params(labelsize=12)
        ax2[i].legend(loc='best', fontsize=15)
    st.pyplot(fig2)

# Espectro
Sa_EO = np.zeros(nregistros*len(T))
Sa_NS = np.zeros(nregistros*len(T))

if open(".\\Espectros.txt",'r').read()=='':    
    #Cálculo de espectros Sa
    for i in range(nregistros):
        Sa_EO[i*len(T):(i+1)*len(T)], T = spectro(Ag_corr_EO[int(i0[i]):int(i0[i+1])],ndatos[i],dt[i],masa=1.0)
        Sa_NS[i*len(T):(i+1)*len(T)], T = spectro(Ag_corr_NS[int(i0[i]):int(i0[i+1])],ndatos[i],dt[i],masa=1.0) 

    np.savetxt('.\\Espectros.txt',np.transpose(np.append(Sa_EO,Sa_NS).reshape((2,len(Sa_EO)))))
else:
    Sa_EO = np.loadtxt('.\\Espectros.txt')[:,0]
    Sa_NS = np.loadtxt('.\\Espectros.txt')[:,1]

#Escalamiento de espectros

Sa = (Sa_EO**2+Sa_NS**2)**0.5
FE = np.zeros(nregistros)

with st.sidebar.container():
    st.write("Factores de Escalamiento")
    for i in range(nregistros):
        FE[i] = st.sidebar.number_input(f'FE {i+1}', 0.1, 50.1, 1.0)

with tab3:
    for i in range(nregistros):
        ax3[i].plot(T,FE[i]*Sa[i*len(T):(i+1)*len(T)], 'blue', linewidth=2,label=(f"registro_FE = {FE[i]:.2f}"))
        ax3[i].plot(T,Sa_e030, 'red', linewidth=2,label='E030')
        ax3[i].set_title(f'Espectro de Pseudo-Aceleraciones {nombre_registros[i]}', fontsize=20, fontweight = 'bold')
        ax3[i].set_xlabel('T (s)', fontsize=15, fontweight = 'bold')
        ax3[i].set_ylabel('Sa (gal)', fontsize=15, fontweight = 'bold')
        ax3[i].set_xlim(0,5)
        ax3[i].tick_params(labelsize=12)
        ax3[i].axvline(0.2*Test, color='m',ls='--', linewidth=0.8,label='0.2T')
        ax3[i].axvline(1.5*Test, color='m',ls='--',linewidth=0.8,label='1.5T')
        ax3[i].legend(loc='best', fontsize=15)
        # ax3[i].set_xscale("log")
    st.pyplot(fig3)


Sa_prom = 0 
for i in range(nregistros):
    Sa_prom += FE[i]*Sa[i*len(T):(i+1)*len(T)]
Sa_prom = Sa_prom/nregistros

with tab4:
    ax4.plot(T,Sa_prom, 'blue', linewidth=2,label="Sa_Promedio")
    ax4.plot(T,Sa_e030, 'red', linewidth=2,label='E030')
    ax4.set_title('Espectro de Pseudo-Aceleraciones', fontsize=20, fontweight = 'bold')
    ax4.set_xlabel('T (s)', fontsize=15, fontweight = 'bold')
    ax4.set_ylabel('Sa (gal)', fontsize=15, fontweight = 'bold')
    ax4.set_xlim(0,5)
    ax4.tick_params(labelsize=12)
    ax4.axvline(0.2*Test, color='m',ls='--', linewidth=0.8,label='0.2T')
    ax4.axvline(1.5*Test, color='m',ls='--',linewidth=0.8,label='1.5T')
    ax4.legend(loc='best', fontsize=15)
    # plt.xscale("log")
    st.pyplot(fig4)

    i = 0
    for archivo in uploaded_files:
        name ='escal_'+archivo.name  
        aux_EO = FE[i]*Ag_corr_EO[int(i0[i]):int(i0[i+1])]
        aux_NS = FE[i]*Ag_corr_NS[int(i0[i]):int(i0[i+1])]
        # np.savetxt('.\\Registros_Escalados\\'+name,
        # np.transpose(np.append(t[int(i0[i]):int(i0[i+1])],np.append(aux_EO,aux_NS)).reshape((3,len(aux_EO)))))
        resultados = np.transpose(np.append(t[int(i0[i]):int(i0[i+1])],np.append(aux_EO,aux_NS)).reshape((3,len(aux_EO))))
        df = pd.DataFrame(resultados,columns=['t (s)','Acc_EO (gal)','Acc_NS (gal)'])

        csv = df.to_csv().encode('utf-8')

        st.download_button(label=f'Descargar {name}',data=csv,file_name='escal_'+archivo.name,mime='text/csv')
        i+=1


# Path('Registros_Escalados').mkdir(exist_ok=True)        
# except:
#     st.write("No se ha cargado ningún archivo")
#     pass
    
    
# ax4.plot(T,Sv, 'blue', linewidth=2,label='registro')
# ax4.plot(T,Sv_e030, 'red', linewidth=2,label='E030')
# ax4.set_title('Espectro de Pseudo-Velocidades', fontsize=20, fontweight = 'bold')
# ax4.set_xlabel('T (s)', fontsize=15, fontweight = 'bold')
# ax4.set_ylabel('Sv (kine)', fontsize=15, fontweight = 'bold')
# ax4.set_xlim(0,5)
# ax4.tick_params(labelsize=12)
# ax4.legend(loc='best', fontsize=15)

# ax5.plot(T,Sd, 'blue', linewidth=2,label='registro')
# ax5.plot(T,Sd_e030, 'red', linewidth=2,label='E030')
# ax5.set_title('Espectro de Desplazamientos', fontsize=20, fontweight = 'bold')
# ax5.set_xlabel('T (s)', fontsize=15, fontweight = 'bold')
# ax5.set_ylabel('D (cm)', fontsize=15, fontweight = 'bold')
# ax5.set_xlim(0,5)
# ax5.tick_params(labelsize=12)
# ax5.legend(loc='best', fontsize=15)



# # Gráficos
# # t, DrX_med, Dr, AgX, Fr, De_pos_f, De_neg_f, Fe_pos_f, Fe_neg_f = plot_graficos(data,masa,fm,ch1,ch2,ch3,ch4,ch5,ch6)
# Data=data
# fm=fm#Frecuencia de Muestreo
# dt = 1/fm
# M=masa*9.81 #Masa en Newton
# t=Data.iloc[:,0]

# AgX=-Data.iloc[:,ch1] #Aceleración absoluta en X en la mesa en g  CH1
# AX=-Data.iloc[:,ch2]  #Aceleración absoluta en X en la estructura en g CH2
# AgY=-Data.iloc[:,ch3] #Aceleración absoluta en Y en la mesa en g  CH3
# #AY=Data[:,4]
# DgX_m=-Data.iloc[:,ch5-1]   #Desplazamiento absoluto en X en la mesa en mm  CH5
# DX_m=-Data.iloc[:,ch6-1]    #Desplazamiento absoluto en X en la estructura en mm CH6
# AY=-Data.iloc[:,ch4-1]      #Aceleración absoluta en Y en la estructura en g CH7
# n=len(t)    #Número de datos

# DrX_med=-DgX_m+DX_m #Desplazamientos relativos MEDIDOS dirección X
# Des_corr=DrX_med[0]
# for j in range(0,n):
#     DrX_med[j] = DrX_med[j] - Des_corr

# #CÁLCULO de los desplazamientos y velocidades por integración numérica
# #Velocidades y Desplzamientos corregidos por línea base
# Vag, Va, Dag, Da = Integ_Corr(AgX,AX,dt,n)

# Fr = AX*M #fuerza cortante cíclica en la estructura
# Dr = Dag-Da #desplazamiento relativo

# #Envolvente de Curva Histerética
# De_neg_f, De_pos_f, Fe_neg_f, Fe_pos_f = Envolvente(DrX_med,Fr,n)

# fig1, ax1 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))
# fig2, ax2 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))
# fig3, ax3 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))


# #return t, DrX_med, Dr, AgX, Fr, De_pos_f, De_neg_f, Fe_pos_f, Fe_neg_f
# tab1, tab2, tab3 = st.tabs(["Structure Relative Displacement", "Ground Acceleration", "Hystheresis"])

# ax1.plot(t,-DrX_med, 'black', linewidth=1, label='Medido')
# ax1.plot(t,-Dr, 'r--', linewidth=1, label='Calculado')
# ax1.set_title('Structure Relative Displacement', fontsize=20, fontweight = 'bold')
# ax1.set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
# ax1.set_ylabel('Displacement (mm)', fontsize=15, fontweight = 'bold')
# ax1.tick_params(labelsize=12)
# ax1.legend(loc='best', fontsize=15)

# with tab1:
#     col1, col2, col3, col4 = st.columns(4)
#     with col1:
#         x = st.number_input('xmín', 0, 100, 0)
#     with col2:
#         X = st.number_input('xmáx', 10, 150, 100)
#     with col3:
#         y = st.number_input('Disp_mín', -10, 0, -5)
#     with col4:
#         Y = st.number_input('Disp_máx', 0, 10, 5)
#         ax1.axis([x, X, y, Y])
#     st.pyplot(fig1)

# ax2.plot(t,AgX, 'black', linewidth=1, label='Base')
# ax2.set_title('Ground Acceleration', fontsize=20, fontweight = 'bold')
# ax2.set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
# ax2.set_ylabel('Acceleration (g)', fontsize=15, fontweight = 'bold')
# ax2.tick_params(labelsize=12)
# ax2.legend(loc='best', fontsize=15)
# # ax1.axis([x, X, y, Y])
# with tab2:
#     col1, col2, col3, col4 = st.columns(4)
#     with col1:
#         x = st.number_input('tmín', 0, 100, 0)
#     with col2:
#         X = st.number_input('tmáx', 10, 150, 100)
#     with col3:
#         y = st.number_input('Ag_mín', -1.5, 0., -0.5)
#     with col4:
#         Y = st.number_input('Ag_máx', 0., 1.5, 0.5)
#         ax2.axis([x, X, y, Y])
#     st.pyplot(fig2)





# ax3.plot(DrX_med, Fr, 'black', linewidth=1, label='Histéresis')
# ax3.plot(De_pos_f,Fe_pos_f, 'r', label='Envolvente')
# ax3.plot(De_neg_f,Fe_neg_f, 'r')
# ax3.set_title('Hystheresis', fontsize=20, fontweight = 'bold')
# ax3.set_xlabel('Displacement X (mm)', fontsize=15, fontweight = 'bold')
# ax3.set_ylabel('Force (N)', fontsize=15, fontweight = 'bold')
# ax3.tick_params(labelsize=12)
# ax3.legend(loc='best', fontsize=15)
# with tab3:
#     # col1, col2, col3, col4 = st.columns(4)
#     # with col1:
#     #     x = st.number_input('xmín', 0, 100, 0)
#     # with col2:
#     #     X = st.number_input('xmáx', 10, 150, 100)
#     # with col3:
#     #     y = st.number_input('ymín', -1, 0, -0.5)
#     # with col4:
#     #     Y = st.number_input('ymáx', 0, 1, 0.5)
#     #     ax3.axis([x, X, y, Y])
#     st.pyplot(fig3)













