#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################################################################
# Asignatura: Sistemas de Comunicación                        ####
# Tema 1: Introducción a la teoría de la señal y los sistemas ####
# Ejercicio 1.1: Simulación de señales y sistemas.            ####
##################################################################
#                                                             ####
# Autor: Fernando Rosa González                               ####
# Fecha: Abril 2015                                           ####
# Area: Teorí­a de la Señal y Comunicaciones                   ####
#                                                             ####
# Las unidades se muestran entre corchetes []                 ####
##################################################################
#                                                             ####
# Descripción: Implementación de la clase de simulación de    ####
#   señales y sistemas.                                          ####
#                                                             ####
##################################################################
import pylab as pyl


class Signal:
    # La clase  Signal permite implementar la simulación de señales
    # estableciendo la ventana de análisis en base a tres parámetros,
    # el instante inicial, el instante final y la frecuencia de muestreo.
    t1 = 0.0                # Instante inicial [s]
    t2 = 0.0                # Instante final   [s]
    fm = 44100.0          # Frecuencia de muestreo [mue/s]
    T = 0.                # Periodo de simulación [s]
    tau = 0.                # Intervalo de muestreo [s]
    N = 0.                # Número total de muestras [mue]
    # Tamaño de la ventana de análisis del espectrograma [mue]
    tam_ven = 1024
    res_frec = 0.         # Resolución en frecuencia [Hz]
    tiempo = pyl.array([])    # Array de tiempos [s]
    frec = pyl.array([])    # Array de frecuencias centradas [Hz]
    frec_pos = pyl.array([])  # Array de frecuencias positivas [Hz]
    ss = pyl.array([])        # Array de la señal [unidades de la señal]
    ff = pyl.array([])     # Array del espectro [unidades de la señal]
    logpl = False            # Estado del plot del espectro

    def __init__(self, t1=0., t2=1., fm=44100.):
        # Los objetos de la clase Signal tienen tres parámetros básicos.
        # t1: instante inicial [s]
        # t2: instante final [s]
        # fm: frecuencia de muestreo [mue/s]
        self.t1 = t1
        self.t2 = t2
        self.fm = fm
        self.T = t2 - t1
        self.tau = 1.0 / self.fm
        self.N = self.fm * self.T
        self.tiempo = pyl.arange(t1, t2, self.tau)
        self.tam_ven = self.N
        self.res_frec = fm / self.tam_ven
        self.frec = pyl.arange(-self.fm / 2., self.fm / 2., self.res_frec)
        self.frec = pyl.concatenate(
            (self.frec[self.N / 2:self.N], self.frec[0:self.N / 2]))
        self.frec_pos = pyl.arange(0., self.fm, self.res_frec)
        self.ss = pyl.zeros(self.N)
        self.ff = pyl.zeros(self.N)
        return(None)

    def step(self, amp_antes=0., amp_desp=1., pos=0.0):
        u'''
        Señal escalon con 'amp_antes' para los tiempos anteriores a pos,
        'amp_desp' para los tiempos posteriores a pos y el valor medio en
        el tiempo 'pos'
        '''
        self.ss = pyl.zeros(self.tiempo.size)
        for n, x in enumerate(self.tiempo):
            if(x > pos):
                self.ss[n] = amp_desp
            elif (x < pos):
                self.ss[n] = amp_antes
            else:
                self.ss[n] = (amp_antes + amp_desp) / 2.
        self.ff = pyl.fft(self.ss) / self.N
        return(None)

    def tone(self, fr=1000., amp=1., fas=0.0):
        u'''
        Señal tono de frecuencia fr, amplitud amp, y fase fas
        '''
        self.ss = amp * pyl.cos(2.0 * pyl.pi * fr * self.tiempo + fas)
        self.ff = pyl.fft(self.ss) / self.N
        return(None)

    def white_noise(self, med=0.0, sig=1.0):
        u'''
        White noise signal with mean med and sigma sig
        '''
        self.ss = pyl.random.normal(med, sig, self.tiempo.size)
        self.ff = pyl.fft(self.ss) / self.N
        return(None)

    def gabor(self, amp=1.0, fr=1000., pos=0.0, sig=1.0, fas=0.0):
        tt = (self.tiempo - pos) / sig
        self.tone(fr, 1.0, fas)
        self.ss *= amp * pyl.exp(-tt * tt)
        self.ff = pyl.fft(self.ss) / self.N
        return(None)

    def pulse(self, amp=1.0, pos=0.5, wid=0.5):
        u'''
        Pulse signal with amplitude, amp, center in time pos, and width, wid.
        '''
        self.ss = pyl.zeros(self.tiempo.size)
        for n, x in enumerate(self.tiempo):
            if((x > pos - wid / 2.) and (x < pos + wid / 2.)):
                self.ss[n] = amp
            elif ((x < pos - wid / 2.) or (x > pos + wid / 2.)):
                self.ss[n] = 0.0
            else:
                self.ss[n] = amp / 2.
        self.ff = pyl.fft(self.ss) / self.N
        return(None)
    u'''
    From here ahead are the operators definitions.
    '''

    def __add__(self, simsig):
        if(type(simsig) == int or type(simsig) ==
           float or type(simsig) == complex):
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss + simsig
        else:
            if (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1):
                print(
                    u"No se pueden sumar señales con diferentes parámetros de simulación")
                return
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss + simsig.ss
        resultado.ff = pyl.fft(resultado.ss) / self.N
        return(resultado)

    __radd__ = __add__

    def __sub__(self, simsig):
        if(type(simsig) == int or type(simsig)
           == float or type(simsig) == complex):
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss - simsig
        else:
            if (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1):
                print(
                    u"No se pueden sumar señales con diferentes parámetros de simulación")
                return
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss - simsig.ss
        resultado.ff = pyl.fft(resultado.ss) / self.N
        return(resultado)

    def __neg__(self):
        resultado = Signal(self.t1, self.t2, self.fm)
        resultado.ss = -self.ss
        resultado.ff = -self.ff
        return(resultado)

    def __mul__(self, simsig):
        if(type(simsig) == int or type(simsig)
           == float or type(simsig) == complex):
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss * simsig
        else:
            if (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1) or \
                    (self.t1 != simsig.t1):
                print(
                    u"No se pueden multiplicar señales con diferentes parámetros de simulación")
                return
            resultado = Signal(self.t1, self.t2, self.fm)
            resultado.ss = self.ss * simsig.ss
        resultado.ff = pyl.fft(resultado.ss) / self.N
        return (resultado)

    __rmul__ = __mul__
    u'''
    From here ahead are the painting methods
    '''

    def pinta(self, titulo=u'sin título', texto=u' '):
        x = self.tiempo
        y = self.ss
        font = {'family': 'serif',
                'color': 'darkred',
                'weight': 'normal',
                'size': 16,
                }
        f, (ax1, ax2) = pyl.subplots(2, 1)
        ax1.plot(x, y, 'b', label=texto)
        ax1.plot([min(x), max(x), pyl.nan, 0.0, 0.0], [
                 0.0, 0.0, pyl.nan, min(y) - 0.1, max(y) + 0.1], 'k--')
        ax1.set_title(titulo, fontdict=font)
        # text(2, 0.65, texto, fontdict=font)
        ax1.set_xlabel('time [s]', fontdict=font)
        ax1.set_ylabel('voltage [V])', fontdict=font)
        ax1.grid(True)
        ax1.set_ylim((min(y) - 0.1, max(y) + 0.1))
        ax1.legend()
        # if(modo=='pos'):
        # off=self.ff;
        # else:
        # off=array(list(self.ff[self.N/2:self.N])+list(self.ff[0:self.N/2]))
        ax2.set_xlabel('frequency [Hz]', fontdict=font)
        if self.logpl is True:
            ax2.plot(self.frec, 20.0 * pyl.log10(abs(self.ff)
                                                 + pyl.spacing(0)))
            ax2.set_ylabel('power [dBV]', fontdict=font)
        else:
            ax2.plot(self.frec, abs(self.ff))
            ax2.set_ylabel('Abs Amplitude [V]', fontdict=font)
        ax2.grid(True)
        pyl.show()

        return(None)

    def specgram(self, tamven=tam_ven):
        self.tam_ven = tamven
        pyl.specgram(self.ss, self.tam_ven, self.fm)
        pyl.show()
        return(None)


class LTI_System:
    fm = 44100.0          # Frecuencia de muestreo [mue/s]
    T = 0.                # Periodo de simulación [s]
    N = 0.                # Número total de muestras [mue]
    res_frec = 0.         # Resolución en frecuencia [Hz]
    tiempo = pyl.array([])    # Array de tiempos [s]
    frec = pyl.array([])    # Array de frecuencias centradas [Hz]
    # Respuesta impulsiva del sistema [Cociente unidades de salida/entrada]
    h = pyl.array([])
    # Función de transferencia [Cociente unidades de salida/entrada]
    H = pyl.array([])

    def __init__(self, T=1.0, fm=44100.0):
        self.fm = fm
        self.T = T
        self.N = T * fm
        self.Tau = 1 / fm
        self.tiempo = pyl.arange(0, self.T, self.Tau)
        self.res_frec = self.fm / self.N
        self.frec = pyl.arange(-self.fm / 2., self.fm / 2., self.res_frec)
        self.frec = pyl.concatenate(
            (self.frec[self.N / 2:self.N], self.frec[0:self.N / 2]))
        self.h = pyl.zeros(self.N)
        self.H = pyl.zeros(self.N, dtype='complex')
        return(None)

    def retraso_temporal(self, t0=0.1):
        u''' Sistema de retraso temporal
        '''
        self.H = pyl.exp(-1.0j * t0 * 2.0 * pyl.pi * self.frec)
        self.h = pyl.ifft(self.H * self.N).real
        return(None)

    def __call__(self, ss):
        out = Signal(t1=ss.t1, t2=ss.t2, fm=ss.fm)
        out.ff = ss.ff * self.H
        out.ss = pyl.ifft(out.ff * out.N).real
        return(out)
