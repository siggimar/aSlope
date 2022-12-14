import matplotlib.pyplot as plt
import matplotlib.colors
from numpy import tanh, arctanh
import scipy.interpolate
import numpy as np

class MODEL_RENDERER():
    def __init__( self, MODEL, STABILITY ) -> None:
        self.colors = {            
            "shear_surface": ( 1, 0, 0 ),
            "background": ( 0.7, 0.7, 0.7 ),
            "contours": ( 0, 0, 0 ),
            "lamellas": ( 0.7, 0.7, 0.7),
            "rock_surface": ( 0.6, 0, 0),
            "gw": (0/255,142/255,194/255), # blue (NPRA)
            "grid": ( 0.2, 0.2, 0.2 ),
            "layers": { # layer rgb, valid for layer_color_max_index < 11
                0: (68/255,79/255,85/255), # dark gray (NPRA)
                1: (93/255,184/255,46/255), # green (NPRA)
                2: (255/255,150/255,0/255), # orange (NPRA)
                3: (237/255,28/255,46/255), # red (NPRA)
                4: (112/255,48/255,160/255), # purple
                5: (0/255,220/255,200/255), # dark cyan
                6: (0/255,120/255,0/255), # dark green
                7: (200/255,120/255,0/255), # dark orange
                8: (180/255,20/255,40/255), # dark red
                9: (80/255,40/255,120/255) # dark purple
            }
        }

        self.lw = {
            "shear_surface": 1,
            "lamellas": 0.1,
            "grid-search": .05,
            "contours": .5
        }

        self.fs_bounds = {
            'low': 1,
            'high': 1.4
        }

        n_contour_interval = (self.fs_bounds['high']-self.fs_bounds['low'])
        n_contour_ints = int( n_contour_interval*10 + 1 )
        n_contour_increment = n_contour_interval / n_contour_ints


        self.contour_ints = [self.fs_bounds['low'] + i*n_contour_increment for i in range(n_contour_ints+1)]

        self.center_map = {}

        self.model = MODEL
        self.fs = STABILITY


    def render( self ):
        figure, ax = plt.subplots( figsize=(10,6) )

        if self.model.GW:
            ax.plot( self.model.GW.x, self.model.GW.y, c=self.colors["gw"] )

        for layer in self.model.layers:
            ax.plot( layer.x, layer.y, c=self.colors["layers"][layer.color_i] )
        ax.plot( self.model.rock.x, self.model.rock.y, c=self.colors["rock_surface"] )


        # plot circular FS
        if self.fs.fs_manager.score_id: # grid search
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list( "", ['red','orange','yellow','green','blue'] )
            
            n_circles = len(self.fs.fs_manager.score_id)

            for key in self.fs.fs_manager.fs: # plot each fs
                fs = self.fs.fs_manager.fs[ key ]
                F = fs['safety']
                c = cmap( self.calc_c_val( F ) ) 
                self.plot_fs( ax, fs, color=c, lw=self.lw['grid-search'], zorder=n_circles-self.fs.fs_manager.score_id[key] )

                # build center map
                cen_key = tuple( list(key[0:2]) )
                if cen_key in self.center_map.keys():
                    if self.center_map[ cen_key ][0] > F:
                        self.center_map[ cen_key ] = [F, c]
                else:
                    self.center_map[ cen_key ] = [F, c]


            # build center map
            x_cen = []
            y_cen = []
            f_cen = []
            for key in self.center_map:
                x, y = list(key)
                x_cen.append(x)
                y_cen.append(y)
                f_cen.append( self.center_map[key][0] )

            n_int = 200
            xi, yi = np.linspace( min(x_cen), max(x_cen), n_int), np.linspace(min(y_cen), max(y_cen), n_int)
            xi, yi = np.meshgrid(xi, yi)

            rbf = scipy.interpolate.Rbf( x_cen, y_cen, f_cen, function='linear' )
            fi = rbf(xi, yi)
            cen_map = ax.pcolormesh( xi, yi, fi, cmap=cmap )
            
            contour_ls = 'dashed'
            c_plt = ax.contour( xi, yi, fi, self.contour_ints, linestyles=[contour_ls], linewidths=[self.lw['contours']], colors=[self.colors['contours']] )
            cb = plt.colorbar( cen_map ) # draw legend
            for some_interval in self.contour_ints: # draw contours on legend
                cb.ax.plot( [0,1], [some_interval]*2, ls=contour_ls, lw=self.lw['contours'], c=self.colors['contours'] )


            key = min(self.fs.fs_manager.score_id, key=self.fs.fs_manager.score_id.get)
            fs = self.fs.fs_manager.fs[ key ]

            if self.fs.fs_manager.grid:
                ax.plot( self.fs.fs_manager.grid['x'], self.fs.fs_manager.grid['y'], c=self.colors["grid"] )

            self.plot_fs( ax, fs, annotate=True, color=(0,0,0) )
            #figure.colorbar(cmap, label='Fs')#, loc='r', length=0.7, ticks=0.1
            #ax.colorbar( cmap )


        else:
            plt_lamellas = True if len(self.fs.fs_manager.fs) == 1 else False
            for key in self.fs.fs_manager.fs:
                fs = self.fs.fs_manager.fs[key]
                self.plot_fs( ax, fs, annotate=True, lamellas=plt_lamellas )

        ax.axis('equal')
        #ax.grid()
        plt.show()


    def calc_c_val( self, F ):
        x_from = arctanh( -0.95 )
        x_to   = arctanh( 0.95 )
        x_scale = x_to - x_from


        f_1 = self.fs_bounds['low']  # ~10%
        f_2 = self.fs_bounds['high'] # ~90%
        f_scale = f_2-f_1
        f_perc = (F-f_1) / f_scale

        x = f_perc * x_scale + x_from

        return tanh( x )/2 + 0.5 # 0-1

    def plot_fs( self, ax, fs, annotate=False, color=None, lw=None, lamellas=False, cmap=None, zorder=9e+20 ):
        
        c = self.colors[ "shear_surface" ]
        if color:
            c = color

        s_lw = self.lw["shear_surface"]
        if lw:
            s_lw = lw


        if annotate:
            # fs_def
            label = 'F'
            if fs['undrained']:
                label += r'$_c$' + r'$_u=$'
            else:
                label += r'$_a$' + r'$_\phi=$'
            if fs['safety']:
                label += str( round(fs['safety'], 3) )

            x_radius = [ fs['x'][0], fs['center_x'], fs['x'][-1] ]
            y_radius = [ fs['y'][0], fs['center_y'], fs['y'][-1] ]
            ax.plot( x_radius, y_radius, c=c, lw=s_lw ) # circle definition
            plt.text( fs['center_x'], fs['center_y']+1, label )#, fontdict=font)


        ax.plot( fs['x'], fs['y'], c=c, lw=s_lw, zorder=zorder ) # shear surface

        if cmap:
            ax.scatter( fs['x'][0], fs['y'][0], c=c, cmap=cmap )


        if lamellas:
            for lam_index in range(len(fs['lamellas'].lamellas)):
                x,y = fs['lamellas'].lamellas[lam_index].get_xy()
                ax.plot( x, y, c=self.colors["lamellas"], lw=self.lw["lamellas"], zorder=-10)