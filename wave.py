# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:32:41 2017

@author: nbutts

This is a simple 2D wave simulator. The generateSurface function sets up the simulation. The 
entire system is visualized with pyqtgraph
"""

import numpy as np
from scipy import signal
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import argparse

def generateSurface(size, A, k):
    '''Returns the simulation surface. In this case it generates a 2D plane. The plane has a 
    tube that goes down the middle and a 90 degree turn at the end. An impulse function
    is generated in the middle of the tube.

    Keyword arguments:
    size -- the size of the simulation surface
    A -- the ampitude of the impulse
    k -- the decay value, lower numbers will cause the wave to continue to propogate
    returns:
    u - the simulation field
    un - the previous simulation field (same as u)
    speed - the propogration speed of the simulation enviornment, essentially the structure
    '''
    un = np.zeros((size+1, size+1))
    u = un
    
    u[size//2-2:size//2+2,0:4] = A
    
    #generate the speed matrix
    speed=np.ones((size+1,size+1)) * k
    
    # Create a pipe
    speed[0:size//2-20,:]=0
    speed[-size//2+20::,0:-5]=0

    # Create a second vent hole
    speed[-size//2+20::,-25:-15]=k
    
    return u, un, speed    

def step_wave_2d(u, un, speed):   
    '''Simulates the 2D system
    Keyword arguments:
    u - the simulation field
    un - the previous simulation field (same as u)
    speed - the propogration speed of the simulation enviornment, essentially the structure
    returns:
    u - the simulation field
    un - the previous simulation field (same as u)
    '''
    size = u.shape[0]
    up=np.zeros((size,size))
    
    uxn = np.vstack((u[1::,:], np.zeros((1, size))))
    uxp = np.vstack((np.zeros((1, size)), u[0:-1, :]))
    uyn = np.hstack((u[:, 1::], np.zeros((size,1))))
    uyp = np.hstack((np.zeros((size, 1)), u[:, 0:-1]))
    
    xcomp = speed * (uxn + uxp - 2 * u)
    ycomp = speed * (uyn + uyp - 2 * u)
    up = 2 * u - un + xcomp + ycomp;
        
    un = u
    u=up

    return u, un


infoStr = """2D Wave simulator"""
parser = argparse.ArgumentParser(description=infoStr)
parser.add_argument('--amp',
                    type=float,
                    default=10,
                    help='The ampitude to use for the initial pulse')
parser.add_argument('--decay',
                    type=float,
                    default=0.1,
                    help='How fast the wave should decay')
parser.add_argument('--dt',
                    type=float,
                    default=0.05,
                    help='How fast the simulation should run in seconds')
parser.add_argument('--size',
                    type=int,
                    default=100,
                    help='How large the simulation surface should be')
parser.add_argument('-s', '--save',
                    action='store_true',
                    help='Save PNG files from the simulation')
args = parser.parse_args()


## Create a GL View widget to display data
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLSurfacePlot')
w.setCameraPosition(distance=150)
w.showMaximized()

## Add a grid to the view
g = gl.GLGridItem()
g.scale(2,2,1)
g.setDepthValue(10)  # draw grid after surfaces since they may be translucent
w.addItem(g)

size = args.size
amp = args.amp
k = args.decay
saveFile = args.save
u, un, speed = generateSurface(size, amp, k)
x = np.linspace(-size/2, size/2, size+1)
p4 = gl.GLSurfacePlotItem(x=x, y = x, shader='heightColor', computeNormals=False, smooth=False)
#p4.shader()['colorMap'] = np.array([0.2, 2, 0.5, 0.2, 1, 1, 0.2, 0, 2])
#p4.translate(10, 10, 0)
w.addItem(p4)

index = 0
def update():
    global p4, z, index, u, un, speed, w, saveFile
    p4.setData(z=u)
    u, un = step_wave_2d(u, un, speed)
    if saveFile:
        saveFile = 'img{0:06}.png'.format(index)
        index = index + 1
        w.grabFrameBuffer().save(saveFile)
    
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(args.dt)

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()


