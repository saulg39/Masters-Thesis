import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.patches as mpatches
from PIL import Image

class BlittedCursor:
    """
    A cross hair cursor using blitting for faster redraw.
    """
    def __init__(self, ax, scale_point):
        self.scale_point = scale_point
        image = Image.open("foobar.png")
        self.width, self.height = image.size
        self.ax = ax
        self.background = None
        self.horizontal_line = ax.axhline(color='k', lw=0.4, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.4, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.72, 0.9, '')
        self._creating_background = False
        self.r=[1,5]
        self.points=[]
        self.pointsx=[]
        self.pointsy=[]
        self.currentline=False
        self.view=[0,0,10]
        self.buildingwidth = 7
        self.redo_grid()
        self.first_clicks=True
        self.first_click=[]
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
            x = round(x/self.r[0],self.r[1])*self.r[0]
            y = round(y/self.r[0],self.r[1])*self.r[0]
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.text.set_position((x+0.1,y+0.1))
            self.text.set_text('(%1.2f, %1.2f)' % (x, y))
            self.ax.figure.canvas.restore_region(self.background)
            if self.currentline!=False:
                cline, = plt.plot([self.currentline[0],x],[self.currentline[1],y],c='0.35',linewidth=1)
                self.ax.draw_artist(cline)
            self.drawing()
            if self.points.count([x,y])!=0:
                self.ax.draw_artist(self.ax.scatter(x,y,c='y',s=20))
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
            x = round(x/self.r[0],self.r[1])*self.r[0]
            y = round(y/self.r[0],self.r[1])*self.r[0]
            if self.first_clicks:
                if self.first_click == []:
                    self.first_click=[x,y]
                else:
                    self.second_click = [x,y]
                    self.first_clicks = False
            else:
                self.points.append([x,y])
                self.pointsx.append(x)
                self.pointsy.append(y)

                
            
            self.ax.figure.canvas.restore_region(self.background)
            self.drawing()
            if self.points.count([x,y])!=0:
                self.ax.draw_artist(self.ax.scatter(x,y,c='y',s=20))
            self.ax.draw_artist(self.horizontal_line)
            self.ax.draw_artist(self.vertical_line)
            self.ax.draw_artist(self.text)
            self.ax.figure.canvas.blit(self.ax.bbox)

            


            self.on_mouse_move("",True)
    def on_key_click(self, event):
        
        
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
                print(self.new_points)
                points = np.array(self.new_points)
                plt.close('all')
                plt.figure(figsize=(11, 8))
                path = evaluate_bezier(points, 50)
                dpath = derivative(path)
                ddpath = derivative(dpath)
                x, y = points[:,0], points[:,1]
                px, py = path[:,0], path[:,1]
                dpx, dpy = dpath[:,0], dpath[:,1]
                ddpx, ddpy = ddpath[:,0], ddpath[:,1]
                plt.plot(px, py, 'b-')
                plt.plot(dpx, dpy, 'k-')
                plt.plot(ddpx, ddpy, 'r-')
                plt.plot(x, y, 'yo')
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



            


scale_point=[21,1]
fig, ax = plt.subplots(figsize=(10,6))
imData = plt.imread("foobar.png")
plt.imshow(imData)
blitted_cursor = BlittedCursor(ax,scale_point)

fig.canvas.mpl_connect('motion_notify_event', blitted_cursor.on_mouse_move)
fig.canvas.mpl_connect('button_press_event', blitted_cursor.on_mouse_click)
fig.canvas.mpl_connect('key_press_event', blitted_cursor.on_key_click)
plt.show()