import streamlit as st
import pandas as pd
# from functions import *
from module import *
import matplotlib.pyplot as plt
from PIL import Image
# plt.style.use('science')


image1 = Image.open('cismid.jpg')
st.image(image1)
st.write('# Escalamiento de Espectros')

image2 = Image.open('uni.jpg')
st.sidebar.image(image2,width=150)
st.sidebar.header('## DATOS')
st.sidebar.markdown("""
[Ejemplo de Archivo de Entrada](https://github.com/elvis1398/Streamlit/blob/main/DA40004.CSV)
""")

# Cargar Archivo
uploaded_file = st.sidebar.file_uploader("Choose a CSV file")
if uploaded_file is not None:
    data = pd.read_csv(uploaded_file,skiprows=14,header=None,engine='python')
else:
    data = pd.read_csv('DA40004.CSV',skiprows=14,header=None,engine='python')
    

# st.write('Datos',data)

# Parámetros
with st.sidebar.container():
    masa = st.number_input('Masa de la Estructura', 20, 70, 50)
    
    fm = st.number_input('Frecuencia de Muestreo', 100, 200, 200)


# Configurar Canales
chanels = (1,2,3,4,5,6,7)
expander1 = st.sidebar.expander("CONFIGURAR CANALES")
ch1 = expander1.selectbox('AgX-Chanel-',chanels,index = 0)
ch2 = expander1.selectbox('AX-Chanel-',chanels,index = 1)
ch3 = expander1.selectbox('AgY-Chanel-',chanels,index = 2)
ch4 = expander1.selectbox('AY-Chanel-',chanels,index = 6)
ch5 = expander1.selectbox('DgX-Chanel-',chanels,index = 4)
ch6 = expander1.selectbox('DX-Chanel-',chanels,index = 5)

# Gráficos
# t, DrX_med, Dr, AgX, Fr, De_pos_f, De_neg_f, Fe_pos_f, Fe_neg_f = plot_graficos(data,masa,fm,ch1,ch2,ch3,ch4,ch5,ch6)
Data=data
fm=fm#Frecuencia de Muestreo
dt = 1/fm
M=masa*9.81 #Masa en Newton
t=Data.iloc[:,0]

AgX=-Data.iloc[:,ch1] #Aceleración absoluta en X en la mesa en g  CH1
AX=-Data.iloc[:,ch2]  #Aceleración absoluta en X en la estructura en g CH2
AgY=-Data.iloc[:,ch3] #Aceleración absoluta en Y en la mesa en g  CH3
#AY=Data[:,4]
DgX_m=-Data.iloc[:,ch5-1]   #Desplazamiento absoluto en X en la mesa en mm  CH5
DX_m=-Data.iloc[:,ch6-1]    #Desplazamiento absoluto en X en la estructura en mm CH6
AY=-Data.iloc[:,ch4-1]      #Aceleración absoluta en Y en la estructura en g CH7
n=len(t)    #Número de datos

DrX_med=-DgX_m+DX_m #Desplazamientos relativos MEDIDOS dirección X
Des_corr=DrX_med[0]
for j in range(0,n):
    DrX_med[j] = DrX_med[j] - Des_corr

#CÁLCULO de los desplazamientos y velocidades por integración numérica
#Velocidades y Desplzamientos corregidos por línea base
Vag, Va, Dag, Da = Integ_Corr(AgX,AX,dt,n)

Fr = AX*M #fuerza cortante cíclica en la estructura
Dr = Dag-Da #desplazamiento relativo

#Envolvente de Curva Histerética
De_neg_f, De_pos_f, Fe_neg_f, Fe_pos_f = Envolvente(DrX_med,Fr,n)

fig1, ax1 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))
fig2, ax2 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))
fig3, ax3 = plt.subplots(constrained_layout=True,facecolor=(1, 1, 1, 1),figsize=(12,8))


#return t, DrX_med, Dr, AgX, Fr, De_pos_f, De_neg_f, Fe_pos_f, Fe_neg_f
tab1, tab2, tab3 = st.tabs(["Structure Relative Displacement", "Ground Acceleration", "Hystheresis"])

ax1.plot(t,-DrX_med, 'black', linewidth=1, label='Medido')
ax1.plot(t,-Dr, 'r--', linewidth=1, label='Calculado')
ax1.set_title('Structure Relative Displacement', fontsize=20, fontweight = 'bold')
ax1.set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
ax1.set_ylabel('Displacement (mm)', fontsize=15, fontweight = 'bold')
ax1.tick_params(labelsize=12)
ax1.legend(loc='best', fontsize=15)

with tab1:
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        x = st.number_input('xmín', 0, 100, 0)
    with col2:
        X = st.number_input('xmáx', 10, 150, 100)
    with col3:
        y = st.number_input('Disp_mín', -10, 0, -5)
    with col4:
        Y = st.number_input('Disp_máx', 0, 10, 5)
        ax1.axis([x, X, y, Y])
    st.pyplot(fig1)

ax2.plot(t,AgX, 'black', linewidth=1, label='Base')
ax2.set_title('Ground Acceleration', fontsize=20, fontweight = 'bold')
ax2.set_xlabel('Time (s)', fontsize=15, fontweight = 'bold')
ax2.set_ylabel('Acceleration (g)', fontsize=15, fontweight = 'bold')
ax2.tick_params(labelsize=12)
ax2.legend(loc='best', fontsize=15)
# ax1.axis([x, X, y, Y])
with tab2:
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        x = st.number_input('tmín', 0, 100, 0)
    with col2:
        X = st.number_input('tmáx', 10, 150, 100)
    with col3:
        y = st.number_input('Ag_mín', -1.5, 0., -0.5)
    with col4:
        Y = st.number_input('Ag_máx', 0., 1.5, 0.5)
        ax2.axis([x, X, y, Y])
    st.pyplot(fig2)





ax3.plot(DrX_med, Fr, 'black', linewidth=1, label='Histéresis')
ax3.plot(De_pos_f,Fe_pos_f, 'r', label='Envolvente')
ax3.plot(De_neg_f,Fe_neg_f, 'r')
ax3.set_title('Hystheresis', fontsize=20, fontweight = 'bold')
ax3.set_xlabel('Displacement X (mm)', fontsize=15, fontweight = 'bold')
ax3.set_ylabel('Force (N)', fontsize=15, fontweight = 'bold')
ax3.tick_params(labelsize=12)
ax3.legend(loc='best', fontsize=15)
with tab3:
    # col1, col2, col3, col4 = st.columns(4)
    # with col1:
    #     x = st.number_input('xmín', 0, 100, 0)
    # with col2:
    #     X = st.number_input('xmáx', 10, 150, 100)
    # with col3:
    #     y = st.number_input('ymín', -1, 0, -0.5)
    # with col4:
    #     Y = st.number_input('ymáx', 0, 1, 0.5)
    #     ax3.axis([x, X, y, Y])
    st.pyplot(fig3)













