import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

colors = ["IndianRed","Green","Purple", "Aqua", "Aquamarine", "Blue", "Bisque", "Black", "BlanchedAlmond"]

def saveResults(fname, populations, values):
    df = pd.DataFrame(values.y.transpose(), columns = populations)
    df.insert(0, 'Time', values.t)
    df.to_csv(fname, float_format='%.5f', sep=',')

def readFile(filename):
    file = open(filename, 'r')
    if (file):        
        return pd.read_csv(filename)
    else: 
        print('Error: it was not possible to open the file!')
        exit

def createFigure(title):
    fig = plt.figure(figsize=(9,7))    
    plt.tick_params()
    plt.title(title)
    plt.xlabel('time (days)')
    plt.ylabel('concentration')
    ax = fig.gca()
    return [fig,ax]
    
def plotFig(id, name, t, u, hasData = False, timeData = None, data = None, logScale= False):
    dir_ = 'figs/'
    fig,ax = createFigure(name)
    if (logScale):
        plt.yscale("log")
    ax.plot(t, u, label=name, color=colors[id%len(colors)])
    if (hasData):
        ax.plot(timeData, data, 'o', color='black')
    ax.legend(loc='best')    
    fig.savefig(dir_ + name + '.svg', format='svg', bbox_inches='tight')
    fig.clear()
    ax.clear()

if __name__ == "__main__":
    odeValues = readFile('ode_output.csv')
    thn = odeValues["Thn"]
    the = odeValues["The"]
    thm = odeValues["Thm"]
    tcd4 = thn + the + thm 
    t = odeValues['t']

    vData = readFile("data/viremia_data.csv")
    plotFig(0,'V',t,odeValues['V'],True,vData["t"],vData["V"])    
    plotFig(1,'Apm',t,odeValues['Apm'])

    il6Data = readFile("data/il6_data.csv")
    plotFig(2,'C',t,odeValues['C'],False,il6Data["t"],il6Data["Il6"])
    
    plotFig(3,'I',t,odeValues['I'])
    plotFig(4,'TCD4',t,tcd4)

    print(odeValues['Tke'])
    plotFig(5,'Tke',t,odeValues['Tke'])

    antibodyData = readFile("data/igm_igg_data.csv")
    plotFig(6,'IgM',t,odeValues['IgM'],False,antibodyData["t"],antibodyData["IgM"])
    plotFig(7,'IgG',t,odeValues['IgG'],False,antibodyData["t"],antibodyData["IgG"])

    plotFig(8,'B',t,odeValues['B'])

    plotFig(9,'Thmi',t,odeValues['Thmi'])

    plotFig(10,'Thn',t,thn)

    plotFig(11,'Thm',t,thm)

    plotFig(12,'The',t,the)