import numpy as np
import matplotlib.pyplot as plt

time = np.array(range(128))/128.0
dt = 1/128.0

def wave(f,phi):
	signal = np.sin(2*np.pi*f*time+phi)
	s1 = np.sin(2*np.pi*1*time+phi)
	s2 = np.sin(2*np.pi*2*time+phi)
	s3 = np.sin(2*np.pi*0*time+phi)
	
	amp = np.fft.fft(signal)/128.0
	powr = amp*np.conjugate(amp)
	
	plt.figure(1)
	plt.subplot(221)
	plt.plot(time,signal)
	plt.plot(time,s1)
	plt.plot(time,s2)
	plt.plot(time,s3)
	plt.scatter(time,signal)
	
	plt.subplot(222)
	freq = np.fft.fftfreq(powr.size,dt)
	freq = np.fft.fftshift(freq)
	#freq = np.roll(freq,63)
	plt.plot(np.log10(powr))
	# freq = [1,5,10,16,32,32.4,64,96,128,132,144,150,192]
	plt.subplot(223)
	plt.plot(np.real(amp))
	plt.subplot(224)
	plt.plot(np.imag(amp))
	plt.show()

def newPlot(f):
	sig = np.sin(2*np.pi*f*time+0.0)
	amp = np.fft.fft(sig)/128.0
	powr = amp*np.conjugate(amp)
	freq = np.fft.fftfreq(powr.size,dt)
	freq = np.fft.fftshift(freq)
	rl = np.real(amp)
	im = np.imag(amp)

	signal = np.sum(im*np.sin(2*np.pi*f*time))+np.sum(rl*np.cos(2*np.pi*f*time))

	#plt.plot(time,signal)
	#plt.show()


















