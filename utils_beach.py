import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Callable
import matplotlib.tri as tri
import matplotlib.animation as animation
from IPython.display import HTML
import matplotlib.cm as cm
import cmocean
import netCDF4 as nc

def getStation(file1, st, comp=1):
    elev1 = []

    with open(file1) as f1:
        l1 = f1.readline()
        l1 = f1.readline()

        # Get timestep info
        info = l1.split()
        skip = float(info[3])
        dt = float(info[2]) / skip

        l1 = f1.readline()

        for lineno, line in enumerate(f1):
            lines = line.split()
            if lines[0] == str(st):
                if comp == 2:
                    x = lines[2]
                else:
                    x = lines[1]

                if x == "NaN":
                    y = np.nan
                else:
                    y = float(x)

                # If dry, set to zero
                if y < -1000:
                    y = 0.

                elev1.append(y)

    if not elev1:
        raise ValueError("Station %d does not exist" % st )

    time = np.arange(len(elev1))*dt*skip/86400.

    return time, np.array(elev1)

def plot_station(st: int,
                data: list,
                title: str='Gauge comparison'):

    plt.clf()
    f = plt.figure()
    
    for file, label in data:
        t, e = getStation(file, st)
        plt.plot(t, e, label=label, linewidth=1)
        
    plt.grid()
    plt.xlabel('Days')
    plt.ylabel('Water level (m)')
    plt.legend()
    plt.title(title)
    plt.show()
    plt.clf()
    return f


def read_triangulation(fname: str) -> (tri.Triangulation, np.ndarray):
    """ Read a fort.14 file and return its triangulation data as a
    matplotlib Triangulation object.
    """
    df = pd.read_csv(fname, delim_whitespace=True, names=list('abcde'),on_bad_lines='skip')

    # Extract the triangulation, and drop them from the main dataframe
    triangles = df.loc[df['e'].notnull()].iloc[:,2:].values
    # Minus one due to Python's indexing
    triangles = triangles.astype(int) - 1
    
    df = df[~df['e'].notnull()]

    # Extract the nodal coordinates
    nodes = df.loc[df['d'].notnull()][['b', 'c']].values.astype(float)

    # Extract the bathymetry
    bathy = df.loc[df['d'].notnull()]['d'].values.astype(float)


    tr = tri.Triangulation(nodes[:,0], nodes[:,1], triangles)
    return tr, bathy


def plot_bathy(fname: str, figsize=None, plot_grid=False, extent=None):
    """ Plot the bathymetry of a specified fort.14 file name. Value is negative below geoid,
    i.e. the values from the fort.14 are inverted.
    """
    plt.clf()
    
    if figsize is None:
        figsize = (10,7)
        
    fig, ax = plt.subplots(figsize=figsize)

    tr, bathy = read_triangulation(fname)
    bathy = -bathy
    ticks = np.linspace(np.min(bathy), np.max(bathy), 50)

    tcf = ax.tricontourf(tr, bathy, levels=ticks)
    if plot_grid:
        ax.triplot(tr)

    if extent is not None:
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])

    fig.colorbar(tcf)
    ax.set_title('Bathymetry (m)')
    plt.show()
    plt.clf()

    return fig
    
def read_63(fname: str, nodes: list, total_nodes: int, col: int=1) -> list:
    """ Given a fort.63-like file, extract the nodal values of the given list of nodes
    and return them as a list. Values from different timestamps are concatenated without
    space, i.e. the last node at timestep n is immediately followed by first node
    at timestep n+1. If the file is a fort.64 and we want the y-velocity values,
    use col=2.
    """
    df = pd.read_csv(fname, skiprows=2, delim_whitespace=True, index_col=False, header=None, names=list('abc'))
    
    # drop the timestamp rows
    ts = list(filter(lambda x: (x % (total_nodes+1) == 0), range(df.shape[0])))
    df = df.drop(ts)
    
    df.iloc[:,0] = df.iloc[:,0].astype(int)

    # Filter the nodes as specified
    cross_elev = df.loc[df.iloc[:,0].isin(nodes)]
    cross_elev.iloc[:,0] = cross_elev.iloc[:,0].astype(int)
    
    return cross_elev.iloc[:,col].values

def read_63_all(fname: str, meshname: str):
    tr,_ = read_triangulation(meshname)
    total_nodes = len(tr.x)
    node_list = list(range(1, total_nodes+1))

    d1 = read_63(fname, nodes=node_list, total_nodes=total_nodes)
    return tr, d1

def make_tri_plot_function(tr: tri.Triangulation, f63: str) -> Callable:
    """ Plot fort.63-like 2D nodal values (e.g. elevation) of an unstructured triangular mesh given
    as a tri.Triangulation object.
    """

    # Use read_63 with all the nodes
    total_nodes = len(tr.x)
    node_list = list(range(1, total_nodes+1))
    z_data = read_63(f63, nodes=node_list, total_nodes=total_nodes)

    fig, ax = plt.subplots(figsize=(10,5))

    def f(frame=0):
        z_data_frame = z_data[total_nodes*frame : total_nodes*(frame+1)]
        tcf = ax.tricontourf(tr, z_data_frame)
        ax.triplot(tr)
        fig.colorbar(tcf)
        #plt.colorbar()
        #plt.colorbar(boundaries=np.linspace(0,1,5)) 

        plt.show()

    return f
    
def read_63_nc(
    filename: str, 
    ) -> np.ndarray:
    """ Load a fort.63 file in netCDF format into a 1D array which contains consecutive nodal values in time.
    Ex: [e1, e2, ..., eN, e1, e2, ..., eN ]
    """
    ds = nc.Dataset(filename)
    
    # Load all timestamps
    data = ds["zeta"][:].data
 
    ds.close()

    return data


def animate_mesh(
    tr: tri.Triangulation,
    data: list,
    plot_grid: bool=False,
    plot_contour: bool=False,
    vmin: float=-1,
    vmax: float=1,
    extent: list=None,
    frames: int=None,
    figsize: tuple=(13,8),
    frameskip: int=1) -> animation:

    """ Given a list of tuples of [ (label, 1D nodal time series data) ], return a function f = f(i) which, for each integer i,
    plots the given data at the ith frame. A 'frame' is simply the values of the nodes
    specified in node_list at a particular timestep. We assume that all the data passed in
    have the same number of nodes and timestep increment.
    """
   
        
    num_subplots = len(data)
    plt.ioff()
    fig, axs = plt.subplots(nrows=1, ncols=len(data), figsize=figsize)
    total_nodes = len(tr.x)
    node_list = list(range(1, total_nodes+1))

    cmap = cmocean.cm.deep.reversed()
    #cmap = cmocean.cm.balance

    cmap.set_bad(color='lightgray')
    #cmap.set_under(color='lightgray')

    tcfs = []
    zs = []
    ticks = np.linspace(vmin, vmax, 70)

    
    for i in range(num_subplots):
        if num_subplots > 1:
            ax = axs[i]
        else:
            ax = axs

        if extent is not None:
            ax.set_xlim([extent[0], extent[1]])
            ax.set_ylim([extent[2], extent[3]])
            
        #z_data = read_63(data[i][1], nodes=node_list, total_nodes=total_nodes)
        z_data = np.copy(data[i][1])

        # Ignore dry nodes
        #z_data[z_data < -1000] = np.nan
        #isbad = z_data < -1000
        #mask = np.any(np.where(isbad[tr.triangles], True, False), axis=1)
        #tr.set_mask(mask)
        
        
        zs.append(z_data)
        z_data_frame = z_data[: total_nodes]
      
            

        tcf = ax.tricontourf(tr, z_data_frame, levels=ticks, extend='both', cmap=cmap, extent=extent, norm='symlog')
        if plot_grid:
            ax.triplot(tr)

      
        cb = fig.colorbar(tcf)
        tcfs.append(tcf)
       

    if frames is None:
        frames = int(len(z_data)/len(node_list))
    
    def animate(frame):
        for i in range(num_subplots):
            if num_subplots > 1:
                ax = axs[i]
            else:
                ax = axs

            ax.clear()
            ax.set_facecolor('lightgray')
            z_data_frame = zs[i][frame*total_nodes: (frame+1)*total_nodes]

           
                
            isbad = z_data_frame < -1000
            mask = np.any(np.where(isbad[tr.triangles], True, False), axis=1)
            tr.set_mask(mask)

            try:
                if plot_contour:
                    ax.tricontour(tr, z_data_frame, levels=10, linewidths=0.5, colors='0.8')
                    
                tcf = ax.tricontourf(tr, z_data_frame, levels=ticks, extend='both', cmap=cmap, extent=extent, norm='symlog')

            except Exception:
                tcf = None
                
            if plot_grid:
                ax.triplot(tr, lw=0.5)

           
            
            ax.set_title(data[i][0])
            
            if extent is not None:
                ax.set_xlim([extent[0], extent[1]])
                ax.set_ylim([extent[2], extent[3]])

            
        
        return tcf

        
    ani = animation.FuncAnimation(fig, animate, frames=range(0,frames,frameskip), interval=200, blit=False)

    return ani


def make_plot_function(
    data: list, 
    node_list: list, 
    xlabel: str='Node number', 
    ylabel: list=['Water elevation (m)'],
    title: str='Sloping beach: y-direction profile') -> Callable:

    """ Given a list of dict of {label: list}, return a function f = f(i) which, for each integer i,
    plots the given data at the ith frame. A 'frame' is simply the values of the nodes
    specified in node_list at a particular timestep. We assume that all the data passed in
    have the same number of nodes and timestep increment.
    """

    def f(frame=0):
        num_subplots = len(data)
        fig, axs = plt.subplots(nrows=1, ncols=len(data), figsize=(13,6))
        M = len(node_list)

        for i in range(num_subplots):
            if num_subplots > 1:
                ax = axs[i]
            else:
                ax = axs
                
            for label, datum in data[i].items():
                line, = ax.plot(node_list, datum[M*frame:M*(frame+1)], '.-', label=label)
    
                # For now, set the limit using the last dataset
                buffer = 0.1
                ax.set_ylim(min(datum), max(datum) + buffer)
        
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel[i])
                ax.set_title(title)
        
                ax.legend()
        plt.show()
        
    return f

def animate63(
    data: list[dict],
    node_list: list,
    xlabel: str='Node number', 
    ylabel: list=['Water elevation (m)'],
    title: str='Sloping beach: y-direction profile',
    frames: int=None) -> animation:

    """ Given a list of dict of {label: datum}, return a function f = f(i) which, for each integer i,
    plots the given data at the ith frame. A 'frame' is simply the values of the nodes
    specified in node_list at a particular timestep. We assume that all the data passed in
    have the same number of nodes and timestep increment.
    """
   
        
    num_subplots = len(data)
    M = len(node_list)
    plt.ioff()
    fig, axs = plt.subplots(nrows=1, ncols=len(data), figsize=(13,6))
    lines = []
    
    for i in range(num_subplots):
        if num_subplots > 1:
            ax = axs[i]
        else:
            ax = axs

        subplot_lines = []
        for label, datum in data[i].items():
            line, = ax.plot(node_list, datum[:M], '.-', label=label)
            subplot_lines.append(line)
    
            ax.legend()
            ax.grid()
            # For now, set the limit using the last dataset
            buffer = 0.1
            ax.set_ylim(min(datum), max(datum) + buffer)
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel[i])
            ax.set_title(title)

        lines.append(subplot_lines)
        

    # lines now contain the following, ex:
    # [ [ADCIRC elev, DG elev], [ADCIRC vel, DG vel] ]
    if frames is None:
        frames = int(len(datum)/len(node_list))
    
    def animate(frame):
        for i in range(num_subplots):
            for (label, datum), line in zip(data[i].items(), lines[i]):
                line.set_data(node_list, datum[M*frame:M*(frame+1)])
        
        return line,
        
    ani = animation.FuncAnimation(fig, animate, frames=range(0,frames,5), interval=100, blit=False)
    return ani

