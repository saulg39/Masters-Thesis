import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import math
import matplotlib.patches as mpatches
from graph_moment import moment_graph
from PIL import Image
import os
import xlrd
from xlutils.copy import copy
from plastic_moment import plastic_moment
from get_data import get_data



class Graphplotter:
    """
    A cross hair cursor using blitting for faster redraw.
    """
    def __init__(self, ax, image, scale_point, shape, row_number, v, k):
        self.scale_point, self.shape, self.v, self.k= [scale_point, shape, v, k]

        self.Material, self.Forming_process, self.Bending_type, self.Actual_Length, self.shape, self.Material_flat, self.Material_corner = get_data(row_number)
        

        self.width, self.height = image.size
        self.ax = ax
        self.background = None
        self.horizontal_line = ax.axhline(color='k', lw=0.4, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.4, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.72, 0.9, '')
        self._creating_background = False
        self.ratio=[1,5]
        self.points=[]
        self.pointsx=[]
        self.pointsy=[]
        self.px=[]
        self.py=[]
        self.redo_grid()
        self.first_clicks=True
        self.first_click = [0,0]
        self.second_click = [0,0]
        self.on_mouse_move("",True)
        



    def redo_grid(self):
        self.ax.set_xlim([0,self.width])
        self.ax.set_ylim([0,self.height])
        ax.set_ylim(ax.get_ylim()[::-1])


    def on_draw(self, event):
        self.create_new_background()

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def create_new_background(self):
        if self._creating_background:
            # discard calls triggered from within this function
            return
        self._creating_background = True
        self.set_cross_hair_visible(False)
        self.ax.figure.canvas.draw()
        self.background = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)
        self.set_cross_hair_visible(True)
        self._creating_background = False

    def drawing(self):
        self.ax.draw_artist(self.ax.scatter([self.first_click[0],self.second_click[0]],[self.first_click[1],self.second_click[1]],c='r',s=30,marker='x'))
        self.ax.draw_artist(self.ax.scatter(self.pointsx,self.pointsy,c='k',s=30,marker='x'))



    def on_mouse_move(self, event, initial=False):
        if self.background is None:
            self.create_new_background()
        if initial ==True:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.restore_region(self.background)
                self.drawing()
                self.ax.figure.canvas.blit(self.ax.bbox)
            return
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.restore_region(self.background)
                self.drawing()
                self.ax.figure.canvas.blit(self.ax.bbox)
        else:
            self.set_cross_hair_visible(True)
            # update the line positions
            x, y = event.xdata, event.ydata
            x = round(x/self.ratio[0],self.ratio[1])*self.ratio[0]
            y = round(y/self.ratio[0],self.ratio[1])*self.ratio[0]
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.text.set_position((x+0.1,y+0.1))
            self.text.set_text('(%1.2f, %1.2f)' % (x, y))
            self.ax.figure.canvas.restore_region(self.background)
            self.drawing()
            if self.points.count([x,y])!=0:
                self.ax.draw_artist(self.ax.scatter(x,y,c='y',s=20))
            if len(self.pointsx)>2.5:
                line, = self.ax.plot(self.px, self.py, 'b-')
                self.ax.draw_artist(line)
            self.ax.draw_artist(self.horizontal_line)
            self.ax.draw_artist(self.vertical_line)
            self.ax.draw_artist(self.text)
            self.ax.figure.canvas.blit(self.ax.bbox)

    def on_mouse_click(self, event):
        if self.background is None:
            self.create_new_background()
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.restore_region(self.background)
                self.drawing()
                self.ax.figure.canvas.blit(self.ax.bbox)
        else:
            x, y = event.xdata, event.ydata
            x = round(x/self.ratio[0],self.ratio[1])*self.ratio[0]
            y = round(y/self.ratio[0],self.ratio[1])*self.ratio[0]
            if self.first_clicks:
                if self.first_click == [0,0]:
                    self.first_click=[x,y]
                else:
                    self.second_click = [x,y]
                    self.first_clicks = False
            else:
                self.points.append([x,y])
                self.pointsx.append(x)
                self.pointsy.append(y)
                if len(self.pointsx)>2.5:
                    self.pointsnp = np.array(self.points)
                    self.path = evaluate_bezier(self.pointsnp, 20)
                    self.px, self.py = self.path[:,0], self.path[:,1]
            
            self.ax.figure.canvas.restore_region(self.background)
            self.drawing()
            if self.points.count([x,y])!=0:
                self.ax.draw_artist(self.ax.scatter(x,y,c='y',s=20))
            if len(self.pointsx)>2.5:
                line, = self.ax.plot(self.px, self.py, 'b-')
                self.ax.draw_artist(line)
            self.ax.draw_artist(self.horizontal_line)
            self.ax.draw_artist(self.vertical_line)
            self.ax.draw_artist(self.text)
            self.ax.figure.canvas.blit(self.ax.bbox)
            self.on_mouse_move("",True)

    def on_mouse_click2(self, event):
        x, y = event.xdata, event.ydata
        print(x,y)
        print(y)
        print(self.Material_flat)


    def on_key_click(self, event):
        if event.key =="c":
            self.points.pop()
            self.pointsx.pop()
            self.pointsy.pop()
        
        if event.key =="o":
            if not self.first_clicks:
                x1 = self.first_click[0]
                y1 = self.first_click[1]
                x0 = self.second_click[0]
                y0 = self.second_click[1]
                self.new_points = []
                self.new_x = []
                self.new_y = []
                for p in self.points:
                    self.new_points.append([self.scale_point[0]*(p[0]-x0)/(x1-x0),self.scale_point[1]*(p[1]-y0)/(y1-y0)])
                self.new_points = np.array(self.new_points)
                plt.close('all')
                fig = plt.figure(figsize=(11, 8))
                self.path = evaluate_bezier(self.new_points, 50)

                self.moment, self.A = moment_graph(shape = self.shape, Material_flat = self.Material_flat, Material_corner= self.Material_corner, last = self.new_points[-1][0])

                plt.plot(self.A, self.moment, 'r-', label='Theoretical')
                dpath = derivative(self.path)
                #ddpath = derivative(dpath)
                x, y = self.new_points[:,0], self.new_points[:,1]
                px, py = self.path[:,0], self.path[:,1]
                dpx, dpy = dpath[:,0], 0.0001*dpath[:,1]
                #ddpx, ddpy = ddpath[:,0], 0.00000001*ddpath[:,1]
                plt.plot(px, py, 'b-', label='Moment')
                plt.ylabel('Moment / KNm')
                plt.xlabel('Curvature / rad')
                plt.grid(True,'both')
                plt.plot(dpx, dpy, 'k-',label='1st Devivative (x10^-4)')
                #plt.plot(ddpx, ddpy, 'r-', label='2nd Devivative (x10^-8)')
                plt.plot(x, y, 'ko')
                fig.canvas.mpl_connect('button_press_event', graph_plotter.on_mouse_click2)
                fig.canvas.mpl_connect('key_press_event', graph_plotter.on_key_click2)
                plt.legend()
                
                plt.show()

    def on_key_click2(self, event):
        plt.close('all')
        fig = plt.figure(figsize=(11, 8))
        if event.key =="u":
            self.Material_flat[5] += 6
            self.moment, self.A = moment_graph(shape = self.shape, Material_flat = self.Material_flat, Material_corner= self.Material_corner, last = self.new_points[-1][0])
        if event.key =="j":
            self.Material_flat[5] -= 6

        if event.key =="i":
            self.Material_flat[5] += 30

        if event.key =="k":
            self.Material_flat[5] -= 30

        if event.key =="q":
            self.Material_flat[2] -= 1

        if event.key =="e":
            self.Material_flat[2] += 1
       
        if event.key =="up":
            self.Material_flat[1] += 6
        if event.key =="down":
            self.Material_flat[1] -= 6
        if event.key =="left":
            self.Material_flat[0] += 2000
        if event.key =="right":
            self.Material_flat[0] -= 2000
        if event.key =="w":
            self.Material_flat[0] += 10000
        if event.key =="x":
            self.Material_flat[1] -= 10000
        if event.key =="a":
            self.A = [i / 1.2 for i in self.A]
        if event.key =="d":
            self.A = [i * 1.2 for i in self.A] 
        self.moment, self.A = moment_graph(shape = self.shape, Material_flat = self.Material_flat, Material_corner= self.Material_corner, last = self.new_points[-1][0])
        plt.plot(self.A, self.moment, 'r-', label='Theoretical')
        #dpath = derivative(self.path)
        #ddpath = derivative(dpath)
        x, y = self.new_points[:,0], self.new_points[:,1]
        px, py = self.path[:,0], self.path[:,1]
        #dpx, dpy = dpath[:,0], dpath[:,1]
        #ddpx, ddpy = ddpath[:,0], ddpath[:,1]
        plt.plot(px, py, 'b-', label='Experimental')
        plt.ylabel('Moment / KNm')
        plt.xlabel('Curvature / rad')
        plt.grid(True,'both')
        #plt.plot(dpx, dpy, 'k-')
        #plt.plot(ddpx, ddpy, 'r-')
        plt.plot(x, y, 'yo')
        fig.canvas.mpl_connect('button_press_event', graph_plotter.on_mouse_click2)
        fig.canvas.mpl_connect('key_press_event', graph_plotter.on_key_click2)
        plt.legend()
        plt.show()

def derivative(points):
    devivative = []
    for i in range(len(points)-1):
        g = (points[i+1][1]-points[i][1])/(points[i+1][0]-points[i][0])
        x = (points[i+1][0]+points[i][0])/2
        devivative.append([x,g])
    return np.array(devivative)

# find the a & b points
def get_bezier_coef(points):
    # since the formulas work given that we have n+1 points
    # then n must be this:
    n = len(points) - 1

    # build coefficents matrix
    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2

    # build points vector
    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]

    # solve system, find a & b
    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2

    return A, B

# returns the general Bezier cubic formula given 4 control points
def get_cubic(a, b, c, d):
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d

# return one cubic curve for each consecutive points
def get_bezier_cubic(points):
    A, B = get_bezier_coef(points)
    return [
        get_cubic(points[i], A[i], B[i], points[i + 1])
        for i in range(len(points) - 1)
    ]

# evalute each cubic curve on the range [0, 1] sliced in n points
def evaluate_bezier(points, n):
    curves = get_bezier_cubic(points)
    return np.array([fun(t) for fun in curves for t in np.linspace(0, 1, n)])



            
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
pic = os.path.join(THIS_FOLDER, 'Graphs','5_normal_strength_6063_T5_aluminium_alloy_beams.png')

fig, ax = plt.subplots(figsize=(15,9))
imData = plt.imread(pic)
image = Image.open(pic)
plt.imshow(imData)
graph_plotter = Graphplotter(ax, image, scale_point=[0.0016,40], shape = "RHS", row_number = (123), v = 0.3, k = -0.46)

fig.canvas.mpl_connect('motion_notify_event', graph_plotter.on_mouse_move)
fig.canvas.mpl_connect('button_press_event', graph_plotter.on_mouse_click)
fig.canvas.mpl_connect('key_press_event', graph_plotter.on_key_click)
plt.show()
