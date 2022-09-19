import numpy as np
from scipy import fftpack, signal
from signal import valid_signals

def Integ_Corr(Ag,A,dt,n):
    '''
    Función que calcula la velocidad y los desplazamientos a partir de la aceleración:
    
    inputs:
    Ag: aceleración del suelo
    A: aceleración absoluta del modelo
    dt: intervalo de tiempo de muestreo
    n: número de datos
    
    outputs:
    Vag: velocidad absoluta y corregida por línea base del suelo
    Va: velocidad absoluta y corregida por línea base del modelo
    Dag: desplazamiento absoluto y corregido por línea base del suelo
    Da: desplazamiento absoluto y corregido por línea base del modelo
    '''
    Vag=np.zeros(n)
    Va=np.zeros(n)
    Dag=np.zeros(n)
    Da=np.zeros(n)
    fm = 1/dt
    
    
    #Cálculo de la velocidad
    for j in range(n-1):
        Vag[j+1] = Vag[j]+(Ag[j]+Ag[j+1])*9810*0.5*dt
        Va[j+1] = Va[j]+(A[j]+A[j+1])*9810*0.5*dt
    #Corrección de la velocidad por linea base
    win=signal.parzen(fm)
    smoothVg=np.zeros(n)
    smoothV=np.zeros(n)
    smoothVg = signal.convolve(Vag, win, mode='same') / sum(win)
    smoothV = signal.convolve(Va, win, mode='same') / sum(win)
    Vag=Vag-smoothVg
    Va=Va-smoothV
    
    #Cálculo del desplazamiento
    for j in range(0,n-1):
        Dag[j+1] = Dag[j]+Vag[j]*dt+(Ag[j]*(dt**2)/3 + Ag[j+1]*(dt**2)/6)*9810
        Da[j+1] = Da[j]+Va[j]*dt+(A[j]*(dt**2)/3 + A[j+1]*(dt**2)/6)*9810
    
    #Corrección de los desplazamientos por linea base
    smoothDg=np.zeros(n)
    smoothD=np.zeros(n)
    smoothDg = signal.convolve(Dag, win, mode='same') / sum(win)
    smoothD = signal.convolve(Da, win, mode='same') / sum(win)
    Dag=(Dag-smoothDg)
    Da=(Da-smoothD)
    
    return Vag, Va, Dag, Da

def Envolvente(Dr_med, Fr, n):
    '''
    Función que calcula la envolvente de una curva histerética:
    
    inputs:
    Dr_med: desplazamiento relativo de la estructura MEDIDO
    Fr: fuerza cortante resultante en la estructura
    n: número de datos
    
    outputs:
    De_neg_f: desplazamiento negativo en la envolvente
    De_pos_f: desplazamiento positivo en la envolvente
    Fe_neg_f: fuerza negativa en la envolvente
    Fe_pos_f: fuerza positiva en la envolvente
    '''
    #ENVOLVENTE
    Fe_pos=np.zeros(n)
    De_pos=np.zeros(n)
    Fe_neg=np.zeros(n)
    De_neg=np.zeros(n)
    k=0
    j=0
    De_max=0
    De_min=0
    m_pos=0.0001 #Valor inicial de la pendiente. Arbitraria pero positiva para que inicie tomando valores
    m_neg=0.0001
    
    #Calculo de la envolvente positiva y negativa
    for i in range(1, n):
        dd=Dr_med[i]-Dr_med[i-1]
        if m_pos>=0:
            if Dr_med[i]>=De_max:
                if dd>0:
                    De_pos[k]=Dr_med[i]
                    Fe_pos[k]=Fr[i]
                    De_max=De_pos[k]
                    k=k+1
        if m_neg>=0:
            if Dr_med[i]<=De_min:
                if dd<0:
                    De_neg[j]=Dr_med[i]
                    Fe_neg[j]=Fr[i]
                    De_min=De_neg[j]
                    j=j+1

    Fe_pos_max=np.amax(Fe_pos)
    Fe_pos_min=np.amin(Fe_neg)

    #Filtrado de envolvente para tomar solo pendientes positivas
    l=1
    De_pos_f=np.zeros(k)
    Fe_pos_f=np.zeros(k)
    De_neg_f=np.zeros(j)
    Fe_neg_f=np.zeros(j)
    m_pos=0
    m_neg=0

    for i in range(1, k):
        if Fe_pos[i]<Fe_pos_max:
            m_pos=(Fe_pos[i]-Fe_pos_f[l-1])/(De_pos[i]-De_pos_f[l-1])
            if m_pos>=0:
                De_pos_f[l]=De_pos[i]
                Fe_pos_f[l]=Fe_pos[i]
                l=l+1
    p=1
    for i in range(1, j):
        #if (De_neg[i]-De_neg[i-1])!=0:
            m_neg=(Fe_neg[i]-Fe_neg_f[p-1])/(De_neg[i]-De_neg_f[p-1])
            if m_neg>=0:
                De_neg_f[p]=De_neg[i]
                Fe_neg_f[p]=Fe_neg[i]
                p=p+1
    #Bucle para obtener arrays que no tengan valores igual a cero al final

    for i in range(l, k):
        De_pos_f[i]=De_pos_f[l-1]
        Fe_pos_f[i]=Fe_pos_f[l-1]

    for i in range(p, j):
        De_neg_f[i]=De_neg_f[p-1]
        Fe_neg_f[i]=Fe_neg_f[p-1]
        
    return De_neg_f, De_pos_f, Fe_neg_f, Fe_pos_f