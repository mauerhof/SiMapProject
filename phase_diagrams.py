# Routines for making phase diagrams
##########################################################################

# rd_f90 - read from f90 generated histogram map
# plotOne - plot one histogram, for one or more snapshots
# plotMany - plot many histograms

##########################################################################
# Plot phase diagram from files (written by phase_diagram.f90)
# Optionally show stacked for many times in Myr (i.e. snapshots)
def plot_phase_diagram(snaps=None, redshifts=None, runDir='./'
                    ,nbins=100, xlog=False, ylog=False
                    ,xvar='nH', yvar='T', wvar='mass'
                    ,filename=None, datDir='./dat',pdf=None
                    ,cmap='viridis',wrg=[1e-6,1e0]
                    ,normalise=False):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    import minirats.utils.py.readwrite.ramses_info as ramses_info
    import minirats.utils.py.readwrite.rd_sim_dir as rd_sim_dir
    import minirats.utils.py.plots as plots
    import matplotlib
    #---------------------------------------------------------------------
    if snaps is None:
        # Find shapshot corresponding to redshifts:
        snaps=rd_sim_dir.getSnaps_redsh(redshifts,dir=runDir)
    plots.initPlot_sCol(title=False)
    fig, ax = plt.subplots()

    plt.minorticks_on()
    xvarInfo=ramses_info.get_varinfo(xvar)
    yvarInfo=ramses_info.get_varinfo(yvar)
    wvarInfo=ramses_info.get_varinfo(wvar)

    hmap=0
    for nout in snaps:                        # Loop over times and stack
        hist=rd_f90(nout, xvar=xvar, yvar=yvar,wvar=wvar
                        ,datDir=datDir,filename=filename)
        hmap=hmap+hist['map']
    # Map now contains the stack over all times and hist contains the last
    # histogram data
    xrg=hist['xrg']
    yrg=hist['yrg']
    #hmap=hmap/hist['wtot']
    #hmap[hmap<1e-10] =1e-10
    if normalise:
        hmap=hmap/np.sum(hmap)

    im=plt.imshow(hmap.T,interpolation='none',aspect='auto' ,origin='lower'
               ,extent=[xrg[0],xrg[1],yrg[0],yrg[1]],cmap=cmap
               ,norm=LogNorm(vmin=wrg[0], vmax=wrg[1]))

    #levels = np.arange(, 1.6, 0.2)

    # Set titles for the axes
    xtitle=xvarInfo['title']
    if hist['xlog']==True: xtitle='log('+xtitle+')'
    plt.xlabel(xtitle)
    ytitle=yvarInfo['title']
    if hist['ylog']: ytitle='log('+ytitle+')'
    plt.ylabel(ytitle)

    xmean=0
    if xmean==1:              # Indicate mean X by weight
        nx=hist['nx']
        ny=hist['ny']
        xval=10.**(np.linspace(xrg[0],xrg[1],nx))
        #yval=10.**(np.linspace(yrg[0],yrg[1],ny))
        nH_map=np.zeros((int(nx),int(ny)))
        for j in range(int(ny)):
           nH_map[:,j-1] = xval[:]
        xw_mean = np.log10(np.sum(hmap*nH_map) / np.sum(hmap))
        plt.plot([xw_mean,xw_mean],[-100,100]
              ,color='firebrick', scalex=False, scaley=False)

    # Show colorbar
    cax = fig.add_axes([0.79, 0.3, 0.02, 0.55])
    cbar=fig.colorbar(im, cax=cax)
    cbar.ax.set_ylabel(wvarInfo['title'])
    if normalise:
        cbar.ax.set_ylabel('P('+wvarInfo['name']+')')
    
    if pdf!=None:
        plots.savePlotPdf(pdf,tight=False)

##########################################################################
# Generate histogram with f90 phase_diagram routine
def generate_f90_phase_diagram(snap=None, redshift=0., RunDir='./'
                    ,nbins=100, xlog=False, ylog=False
                    ,xrg=[0,1], yrg=[0,1]
                    ,xvar='nH', yvar='T', wvar='mass'
                    ,filename=None, datDir='./dat'):
    import os
    import sys
    #from pathlib import Path
    from minirats.utils.py.readwrite import rd_sim_dir
    from minirats.utils.py.readwrite import ramses_info
    #------------------------------------------------------------------------
    if snap is None:
        # Find shapshot corresponding to redshift:
        snap = rd_sim_dir.getSnaps_redsh(redshift,dir=RunDir)[0]
    snapStr="%5.5i"%(snap)
    print('Reading snapshot ',snap)
    ixlog=0
    iylog=0
    irdRT=0
    if xlog is True:
        ixlog=1
    if ylog is True:
        iylog=1
    cwd = os.getcwd()
    if filename is None:
        filename= "%s/%s/pH_%s_%s_%s_%s.dat" \
                                      %(cwd,datDir,xvar,yvar,wvar,snapStr)
    cmd = '$MINIRATS/utils/f90/standalone/phase_diagram'                 \
      +' -inp '+snapStr                                                  \
      +' -xvr '+xvar+' -yvr '+yvar                                       \
      +' -xlg '+str(ixlog)+' -ylg '+str(iylog)                           \
      +' -vxd '+str(xrg[0])+' -vxu '+str(xrg[1])                         \
      +' -vyd '+str(yrg[0])+' -vyu '+str(yrg[1])                         \
      +' -wvr '+wvar                                                     \
      +' -nbi '+str(nbins)                                               \
      +' -out '+filename                                                 \
      +' -rt  '+str(irdRT)

    print(cmd)
    os.chdir(RunDir)
    os.system(cmd)
    os.chdir(cwd)

    

      #  +' -inp output_'+snapStr                                            \

##########################################################################
# Read and return histogram from f77-generated file
def rd_f90(nout, xvar='nH', yvar='T', wvar='mass'
                      ,filename=None,datDir='./dat',outsideRegion=False):
    import numpy as np
    #------------------------------------------------------------------------
    snapStr="%5.5i"%(nout)
    if filename is None:
        filename= "%s/pH_%s_%s_%s_%s.dat" \
                                      %(datDir,xvar,yvar,wvar,snapStr)
    f=open(filename,'rb')

    recl = np.zeros(1,dtype=np.int32)
    tmp = np.fromfile(f, dtype='int32', count=1)
    nx = int(np.fromfile(f, dtype='int32', count=1))
    ny = int(np.fromfile(f, dtype='int32', count=1))
    tmp = np.fromfile(f, dtype='int32', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    print('nx,ny=',nx,ny)

    ixlog = np.fromfile(f, dtype='int32', count=1)
    xrg0 = np.fromfile(f, dtype='float64', count=1)
    xrg1 = np.fromfile(f, dtype='float64', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    #print('ixlog,xrg0,xrg1=',ixlog,xrg0,xrg1)
    
    iylog = np.fromfile(f, dtype='int32', count=1)
    yrg0 = np.fromfile(f, dtype='float64', count=1)
    yrg1 = np.fromfile(f, dtype='float64', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    #print('iylog,yrg0,yrg1=',iylog,yrg0,yrg1)
    
    wtot = np.fromfile(f, dtype='float64', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    tmp = np.fromfile(f, dtype='int32', count=1)
    #print('wtot=',wtot)

    maptmp = np.fromfile(f, dtype='float64', count=nx*ny)
    #print('minmax=',len(maptmp),min(maptmp),max(maptmp))
    map = np.reshape(maptmp,(nx,ny),order='F')
    #print('minmax=',np.shape(map),np.min(map),np.max(map))
    f.close()

    xrg=[xrg0[0],xrg1[0]]
    yrg=[yrg0[0],yrg1[0]]
    xlog=False
    if ixlog==1: xlog=True
    ylog=False
    if iylog==1: ylog=True

    return {'map':map, 'nx':nx, 'ny':ny, 'xrg':xrg, 'yrg':yrg
            ,'wtot':wtot,'xlog':xlog, 'ylog':ylog}
            


if __name__ == '__main__':
	
	generate_f90_phase_diagram(12)



