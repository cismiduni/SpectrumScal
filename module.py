import numpy as np
from scipy import fftpack, signal
from signal import valid_signals

def intsignal(dt0, Accg0, ndatos0):
    ''' Función que interpola el registro de aceleraciones para una mayor discretización de datos
    input
    dt0:    intervalo de tiempo del registro
    Accg0:  registro de aceleraciones

    retorna
    Accg:   registro de aceleraciones interpolado
    ndatos: número de puntos de Accg
    dt:     nuevo intervalo de tiempo
    '''
    #Preproces. señal aumentar puntos dt (fmuestreo Hz)
    # Nyquist
    dt = dt0/4
    nint = int(dt0/dt) # número de intervalos de interpolación
    ndatos = (ndatos0-1)*nint + 1
    Accg = np.zeros(ndatos)

    Accg[ndatos-1] = Accg0[ndatos0-1]
    
    k = 0
    for i in range(ndatos0-1):
        for j in range(nint):
            Accg[k] = (Accg0[i]*(nint-j) + Accg0[i+1]*j)/nint
            k += 1
    return Accg, ndatos, dt

def NewmarkBeta(ndatos, Accg, dt, masa, T, h):
    ''' Función que aplica el método de Newmark beta a un registro de aceleraciones
    input
    Accg:       registro de aceleraciones
    ndatos:     tamaño del registro
    dt:         intervalo de tiempo del registro
    masa:       masa de la estructura
    T:          periodo de la estructura
    h:          factor de amortiguamiento de la estructura

    retorna
    amax: Aceleración absoluta máxima
    '''
    
    w = 2*np.pi/T #frecuencia angular
    c = 2*h*w*masa
    k = masa*w**2 #rigidez

    #valores iniciales de la respuesta
    ac, ve, de = [-Accg[0],0.0,0.0]
    aAcc = np.zeros(ndatos) #absolute acceleration
    
    beta = 1/6
    coef1, coef2, coef3, coef4, coef5 = [1/(2*beta*dt),1/(beta*dt**2),1/(beta*dt),1/(2*beta),1/(4*beta)-1]

    #Newmark's Beta Method
    for i in range(1, ndatos):
        dag = Accg[i] - Accg[i-1] #variación de aceleraciones del suelo
        dk = k + coef1*c + coef2*masa #variación de rigidez
        dfuerza = -masa*dag + masa*(coef3*ve + coef4*ac) + c*(coef4*ve + coef5*ac*dt)
        
        dd = dfuerza/dk #variación del desplazamiento
        dv = coef1*dd - coef4*ve - coef5*ac*dt #variación de la velocidad
        da = coef2*dd - coef3*ve - coef4*ac #variación de la aceleración
        
        ac, ve, de = [da + ac, dv + ve, dd + de]
        aAcc[i] = ac + Accg[i]

    
    amax = np.amax(np.absolute(aAcc))
    return amax

def Filt_Corr(Ag,dt,n,f1=0.1,f2=25):
    '''
    Función que calcula la velocidad y los desplazamientos a partir de la aceleración:
    
    inputs:
    Ag: aceleración del suelo
    dt: intervalo de tiempo de muestreo
    n: número de datos
    f1: 0.1 #frecuencia mínima en Hz
    f2: 25 #frecuencia máxima en Hz
    
    outputs:
    Ag: Aceleración corregida
    '''
    
    fm = 1/dt
    #Filtro de señal por paso de banda
    sos = signal.butter(4, [f1,f2], btype='bandpass', fs=fm, output='sos') 
    Ag = signal.sosfilt(sos, Ag)
    
    
    #Corrección de la velocidad por linea base
    win=signal.parzen(fm)
    smoothAg=np.zeros(n)
    smoothAg = signal.convolve(Ag, win, mode='same') / sum(win)
    Ag=Ag-smoothAg
        
   
    return Ag


def spectro(Accg,ndatos,dt,masa):
    T = np.arange(0.0,5.01,0.01,dtype = 'double')
    T[0] = 0.02
    nT = T.shape[0]

    h = 0.05

    Sa = np.zeros(nT)
    Sv = np.zeros_like(Sa)
    Sd = np.zeros_like(Sa)

    for i2 in range(0, nT):
        periodo = T[i2]
        amax = NewmarkBeta(ndatos, Accg, dt, masa, periodo, h)
        Sa[i2] = amax
        # Sv[i2,i1] = vmax
        # Sd[i2,i1] = dmax
            
    return Sa, T

