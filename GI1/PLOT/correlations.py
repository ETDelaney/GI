import matplotlib.pyplot as plt
import numpy as np
import parameters

def correlations(rec0=0,rec1=1,effective=True):
	"""
	correlations(idx1=0,idx2=1,effective=True)

	Plot unprocessed, processed and effective correlation functions.

	INPUT:
	------
	rec0, rec1:		indeces of the receivers used in the correlation. 
	effective:		plot effective correlation. Must already be stored in OUTPUT/correlations	

	OUTPUT:
	-------
	none
	
	Last updated: 24 March 2016.
	"""

	#==============================================================================
	#- Input.
	#==============================================================================

	#- Load frequency and time axes. ----------------------------------------------

	fn='OUTPUT/correlations/f'
	fid=open(fn,'r')
	f=np.load(fid)
	fid.close()

	f[0]=0.1*f[1]

	df=f[1]-f[0]
	n=len(f)

	fn='OUTPUT/correlations/t'
	fid=open(fn,'r')
	t=np.load(fid)
	fid.close()

	dt=t[1]-t[0]

	#- Load frequency-domain correlation functions. -------------------------------

	fn='OUTPUT/correlations/ccf_'+str(rec0)+'_'+str(rec1)
	fid=open(fn,'r')
	ccf=np.load(fid)
	fid.close()

	fn='OUTPUT/correlations/ccf_proc_'+str(rec0)+'_'+str(rec1)
	fid=open(fn,'r')
	ccf_proc=np.load(fid)
	fid.close()

	if effective==True:
		fn='OUTPUT/correlations/ccf_eff_'+str(rec0)+'_'+str(rec1)
		fid=open(fn,'r')
		ccf_eff=np.load(fid)
		fid.close()

	#- Read receiver positions and compute frequency-domain Green function. -------

	p=parameters.Parameters()

	d=np.sqrt((p.x[rec0]-p.x[rec1])**2 + (p.y[rec0]-p.y[rec1])**2)

	#==============================================================================
	#- Time-domain correlations.
	#==============================================================================

	dummy=np.real(np.fft.ifft(ccf)/dt)
	cct=np.zeros(np.shape(dummy))
	cct[0.5*n:n]=dummy[0:0.5*n]
	cct[0:0.5*n]=dummy[0.5*n:n]

	dummy=np.real(np.fft.ifft(ccf_proc)/dt)
	cct_proc=np.zeros(np.shape(dummy))
	cct_proc[0.5*n:n]=dummy[0:0.5*n]
	cct_proc[0:0.5*n]=dummy[0.5*n:n]

	if effective==True:
		dummy=np.real(np.fft.ifft(ccf_eff)/dt)
		cct_eff=np.zeros(np.shape(dummy))
		cct_eff[0.5*n:n]=dummy[0:0.5*n]
		cct_eff[0:0.5*n]=dummy[0.5*n:n]

	#==============================================================================
	#- Plot results.
	#==============================================================================

	scale_proc=np.max((np.max(np.abs(cct_proc)), np.max(np.abs(cct_eff))))
	scale=np.max(np.abs(cct))

	if (1==1):

		plt.plot(t,cct_proc/scale_proc,'r',linewidth=1.5)
		if effective==True: plt.plot(t,cct_eff/scale_proc,'b',linewidth=1.5)
		plt.plot(t,cct/scale,'--k',linewidth=1.5)
	
		plt.text(0.5*np.max(t),0.7,'x %0.2g' % scale,color=(0.1,0.1,0.1),fontsize=14)
		plt.text(0.5*np.max(t),0.6,'x %0.2g' % scale_proc,color=(0.9,0.1,0.1),fontsize=14)
		plt.text(0.5*np.max(t),0.5,'x %0.2g' % scale_proc,color=(0.1,0.1,0.9),fontsize=14)

		plt.xlim((0.7*np.min(t),0.7*np.max(t)))
		plt.ylim((-1.1,1.1))

		plt.title('Correlation functions (black=original, red=processed, blue=effective)')
		plt.xlabel('time [s]')
		plt.ylabel('Green functions [s/kg], [s/k]*unit(T)')
		plt.show()

	if (1==0):

		if effective==True: 
			plt.plot(t,(cct_proc-cct_eff)/scale_proc,'--',color=(0.8,0.8,0.8),linewidth=1.5)
			plt.plot(t,cct_eff/scale_proc,'r',linewidth=1.5)
			
		plt.plot(t,cct_proc/scale_proc,'k',linewidth=2.0)
	
		plt.text(0.25*np.max(t),-0.9,'x %0.2g' % scale_proc,color=(0.1,0.1,0.1),fontsize=14)

		plt.xlim((0.35*np.min(t),0.35*np.max(t)))
		plt.ylim((-1.1,1.1))

		plt.title('Correlation functions (black=processed, red=effective, dashed=error)')
		plt.xlabel('time [s]')
		plt.ylabel('Green functions [s/kg], [s/k]*unit(T)')
		plt.show()
