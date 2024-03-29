import numpy as np
import matplotlib.pyplot as plt
import source as s
import green as g
import parameters
import time

def snapshot(rec=0,t=1.0, mg_level=5, mg_tol=0.05, minvalplot=0.0, maxvalplot=0.0, plot=0, save=0, verbose=0, precomputed=0, dir_precomputed='OUTPUT/'):

	"""
	
	snapshot(rec=0,t=1.0, mg_level=5, mg_tol=0.05, minvalplot=0.0, maxvalplot=0.0, plot=0, save=0, verbose=0, precomputed=0, dir_precomputed='OUTPUT/')

	Compute and plot correlation wavefield.

	INPUT:
	------

	rec:				index of receiver.
	t: 					time in s.
	mg_level:			level for multi-grid solver.
	mg_tol:				tolerance for multi-grid solver.
	minvalplot:			minimum of colour scale, ignored when 0.
	maxvalplot: 		maximum of colour scale, ignored when 0.
	plot:				plot when 1.
	save:				save as png when 1.
	verbose:			give screen output when 1.
	precomputed:		set to 1 if precomputed frequency-domain correlation field available (see precompute() below).
	dir_precomputed:	directory where precomputed correlation field is located.

	OUTPUT:
	-------

	C:		2D time-domain correlation wavefield [N^2 s / m^4].
	x,y:	2D axes [m].

	Last updated: 27 May 2016.
	"""

	#==============================================================================
	#- Initialisation.
	#==============================================================================

	p=parameters.Parameters()

	#- Spatial grid.
	x_line=np.arange(p.xmin,p.xmax,p.dx)
	y_line=np.arange(p.ymin,p.ymax,p.dy)
	x,y=np.meshgrid(x_line,y_line)

	#- Frequency line.
	f=np.arange(p.fmin-p.fwidth,p.fmax+p.fwidth,p.df)
	omega=2.0*np.pi*f

	#- Power-spectral density.
	S,indeces=s.space_distribution()
	instrument,natural=s.frequency_distribution(f)
	filt=natural*instrument*instrument

	C=np.zeros(np.shape(x))

	#==============================================================================
	#- Load forward interferometric wavefields.
	#==============================================================================

	if precomputed==1:

		fn=dir_precomputed+'/cf_'+str(rec)
		fid=open(fn,'r')
		Cfull=np.load(fid)
		fid.close()

		for idx in range(len(x_line)):
			for idy in range(len(y_line)):

				C[idy,idx]=np.real(np.sum(Cfull[idy,idx,:]*np.exp(1j*omega*t)))*p.df


	#==============================================================================
	#- Compute correlation field for a specific time.
	#==============================================================================

	else:

		#- First multi-grid stage. ------------------------------------------------

		if (verbose==1): print 'First multi-grid stage'

		#- March through the spatial grid.
		for idx in range(0,len(x_line),mg_level):

			if (verbose==1): print str(100*float(idx)/float(len(x_line)))+' %'

			for idy in range(0,len(y_line),mg_level):

				C_proto=np.zeros(len(omega),dtype=complex)

				#- March through all sources.
				for k in indeces:

					C_proto+=g.conjG1_times_G2(x[idy,idx],y[idy,idx],p.x[rec],p.y[rec],x[k],y[k],omega,p.dx,p.dy,p.rho,p.v,p.Q)*S[k]
					C_proto=np.conj(C_proto)

				#- Transform to time domain.
				C[idy,idx]=np.real(np.sum(filt*C_proto*np.exp(1j*omega*t)))

		#- Second multi-grid stage. -----------------------------------------------

		if (verbose==1): print 'Second multi-grid stage'

		c_max=np.max(np.abs(C))

		#- March through the spatial grid.
		for idx in range(mg_level,len(x_line)-mg_level):

			if (verbose==1): print str(100*float(idx)/float(len(x_line)))+' %'

			for idy in range(mg_level,len(y_line)-mg_level):

				if (np.max(np.abs(C[(idy-mg_level):(idy+mg_level),(idx-mg_level):(idx+mg_level)]))>mg_tol*c_max):

					C_proto=np.zeros(len(omega),dtype=complex)

					#- March through all sources.
					for k in indeces:

						C_proto+=g.conjG1_times_G2(x[idy,idx],y[idy,idx],p.x[rec],p.y[rec],x[k],y[k],omega,p.dx,p.dy,p.rho,p.v,p.Q)*S[k]
						C_proto=np.conj(C_proto)

					#- Transform to time domain.
					C[idy,idx]=np.real(np.sum(filt*C_proto*np.exp(1j*omega*t)))

		#- Normalisation.
		C=C*p.dx*p.dy*p.df

	#==============================================================================
	#- Plot.
	#==============================================================================

	if (plot==1 or save==1):

		if (minvalplot==0.0 and maxvalplot==0.0):
			maxvalplot=0.8*np.max(np.abs(C))
			minvalplot=-maxvalplot

		font = {'size'   : 14,}

		#- Plot interferometric wavefield. ----------------------------------------

		plt.pcolor(x,y,C,cmap='RdBu',vmin=minvalplot,vmax=maxvalplot)

		#- Plot receiver positions. -----------------------------------------------

		for k in range(p.Nreceivers):

			plt.plot(p.x[k],p.y[k],'kx')
			plt.text(p.x[k]+3.0*p.dx,p.y[k]+3.0*p.dx,str(k),fontdict=font)

		plt.plot(p.x[rec],p.y[rec],'ro')
		plt.text(p.x[rec]+3.0*p.dx,p.y[rec]+3.0*p.dx,str(rec),fontdict=font)

		#- Embellish the plot. ----------------------------------------------------

		plt.colorbar()
		plt.axis('image')
		plt.title('correlation field, t='+str(t)+' s')
		plt.xlim((p.xmin,p.xmax))
		plt.ylim((p.ymin,p.ymax))
		plt.xlabel('x [m]')
		plt.ylabel('y [m]')

		if (plot==1):
			plt.show()
		if (save==1):
			fn='OUTPUT/'+str(t)+'.png'
			plt.savefig(fn)
			plt.clf()

	#==============================================================================
	#- Return.
	#==============================================================================
	
	return C, x, y


def movie(time_axis, minvalplot=0.0, maxvalplot=0.0, verbose=0):

	"""
	movie(time_axis, mg_level=5, mg_tol=0.05, minvalplot=0.0, maxvalplot=0.0, verbose=0)

	Compute correlation wavefield and save png figures to /OUTPUT.

	INPUT:
	------

	time_axis: 	array containing time values in s for which figures will be saved.
	minvalplot:	minimum of colour scale, ignored when 0.
	maxvalplot: maximum of colour scale, ignored when 0.
	verbose:	give screen output when 1.

	OUTPUT:
	-------

	Snapshots of the interferometric wavefield saved as png files to /OUTPUT

	Last updated: 5 May 2016.
	"""

	#==============================================================================
	#- Initialisation.
	#==============================================================================

	p=parameters.Parameters()

	#- Spatial grid.
	x_line=np.arange(p.xmin,p.xmax,p.dx)
	y_line=np.arange(p.ymin,p.ymax,p.dy)
	x,y=np.meshgrid(x_line,y_line)

	#- Frequency line.
	f=np.arange(p.fmin-p.fwidth,p.fmax+p.fwidth,p.df)
	omega=2.0*np.pi*f

	#- Power-spectral density.
	S,indeces=s.space_distribution()
	instrument,natural=s.frequency_distribution(f)
	filt=natural*instrument*instrument

	#- Check if the indeces are actually available. If not, interrupt.
	if len(indeces)==0:
		print 'Correlation field cannot be computed because source index array is empty.'
		return

	#==============================================================================
	#- Compute correlation field for specific times and store.
	#==============================================================================

	C=np.zeros((len(y_line),len(x_line),len(time_axis)))

	#- March through the spatial grid. --------------------------------------------
	for idx in range(len(x_line)):

		if (verbose==1):
			print str(100*float(idx)/float(len(x_line)))+' %'

		for idy in range(len(y_line)):

			C_proto=np.zeros(len(omega),dtype=complex)

			#- March through all sources.
			for k in indeces:

				C_proto+=g.conjG1_times_G2(x[idy,idx],y[idy,idx],p.x[1],p.y[1],x[k],y[k],omega,p.dx,p.dy,p.rho,p.v,p.Q)*S[k]
				C_proto=np.conj(C_proto)

			#- Transform to time domain.

			for t in range(len(time_axis)):

				C[idy,idx,t]=np.real(np.sum(filt*C_proto*np.exp(1j*omega*time_axis[t])))


	#- Normalisation.
	C=C*p.dx*p.dy*p.df

	#==============================================================================
	#- Save images.
	#==============================================================================

	for t in range(len(time_axis)):

		if (minvalplot==0.0 and maxvalplot==0.0):
			maxvalplot=0.8*np.max(np.abs(C[:,:,t]))
			minvalplot=-maxvalplot

		plt.pcolor(x,y,C[:,:,t],cmap='RdBu',vmin=minvalplot,vmax=maxvalplot)

		font = {'family' : 'sansserif', 'color'  : 'darkred', 'weight' : 'normal', 'size'   : 14,}
		plt.plot(p.x[0],p.y[0],'ko')
		plt.plot(p.x[1],p.y[1],'ko')
		plt.text(0.9*p.x[0],p.y[0],'1',fontdict=font)
		plt.text(1.1*p.x[1],p.y[1],'2',fontdict=font)
		
		plt.colorbar()
		plt.axis('image')
		plt.title('correlation field, t='+str(time_axis[t])+' s')
		plt.xlabel('x [m]')
		plt.ylabel('y [m]')

		fn='OUTPUT/'+str(time_axis[t])+'.png'
		plt.savefig(fn)
		plt.clf()


def precompute(rec=0,verbose=0):

	"""
	precompute(verbose=0)

	Compute correlation wavefield in the frequency domain and store for in /OUTPUT for re-use in kernel computation.

	INPUT:
	------

	rec:		index of reference receiver.
	verbose:	give screen output when 1.

	OUTPUT:
	-------

	Frequency-domain interferometric wavefield stored in /OUTPUT.

	Last updated: 16 May 2016.
	"""

	#==============================================================================
	#- Initialisation.
	#==============================================================================

	p=parameters.Parameters()

	#- Spatial grid.
	x_line=np.arange(p.xmin,p.xmax,p.dx)
	y_line=np.arange(p.ymin,p.ymax,p.dy)
	x,y=np.meshgrid(x_line,y_line)

	#- Frequency line.
	f=np.arange(p.fmin-p.fwidth,p.fmax+p.fwidth,p.df)
	omega=2.0*np.pi*f

	#- Power-spectral density.
	S,indeces=s.space_distribution()
	instrument,natural=s.frequency_distribution(f)
	filt=natural*instrument*instrument

	#- Check if the indeces are actually available. If not, interrupt.
	if len(indeces)==0:
		print 'Correlation field cannot be computed because source index array is empty.'
		return

	C=np.zeros((len(y_line),len(x_line),len(omega)),dtype=complex)

	#==============================================================================
	#- Compute correlation field for all frequencies and store.
	#==============================================================================

	#- March through the spatial grid. --------------------------------------------
	for idx in range(len(x_line)):

		if (verbose==1): print str(100*float(idx)/float(len(x_line)))+' %'
		t1=time.time()

		for idy in range(len(y_line)):

			#- March through all sources.
			for k in indeces:

				C[idy,idx,:]+=S[k]*filt*g.conjG1_times_G2(x[idy,idx],y[idy,idx],p.x[rec],p.y[rec],x[k],y[k],omega,p.dx,p.dy,p.rho,p.v,p.Q)
				C[idy,idx,:]=np.conj(C[idy,idx,:])

		t2=time.time()
		print t2-t1

	#- Normalisation.
	C=C*p.dx*p.dy

	#==============================================================================
	#- Save interferometric wavefield.
	#==============================================================================

	fn='OUTPUT/cf_'+str(rec)
	fid=open(fn,'w')
	np.save(fid,C)
	fid.close()