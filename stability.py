from tqdm import tqdm
from lamelle import lamellas as lam
from numpy import cos, arccos, arctan2
import copy


class sircular_cylindric_fs():
    def __init__( self, model, n_lamelle=30 ):
        self.model = model
        self.n_lamelle = n_lamelle
        self.fs_manager = sircular_cylindric_fs.failure_surface_manager()

    def sort( self ):
        self.fs_manager.sort()

    def print_geom( self ):
        for fs_key in self.fs_manager.fs:
            fs = self.fs_manager.fs[fs_key]
            lamellas = fs["lamellas"].lamellas
            for lamella in lamellas:
                coords = ''
                coords += str(lamella.x1_4) + ", " + str(lamella.y1)
                coords += ", " + str(lamella.x2_3) + ", " + str(lamella.y2)
                coords += ", " + str(lamella.x2_3) + ", " + str(lamella.y3)
                coords += ", " + str(lamella.x1_4) + ", " + str(lamella.y4)
                print( coords )


    def print_gv( self ): # for backcalculations
        if self.model.GW:
            for fs_key in self.fs_manager.fs:
                fs = self.fs_manager.fs[fs_key]
                lamellas = fs["lamellas"].lamellas
                i=0
                for lamella in lamellas:
                    i += 1
                    out = str( i ) + ': '

                    h_w = lamella.u / self.model.GW.gamma
                    h = max( 0, lamella.y_cen + h_w )
                    out += str( h )
                    print( out )
        else:
            print('no GW defined')


    def grid_search( self, search_field=None, increments=None, clear_FS=False ):
        if clear_FS:
            self.fs_manager.clear_fs()

        # set grid
        if search_field: # implement search box
            x_from = search_field[ 'x_from' ]
            y_from = search_field[ 'y_from' ]
            x_to = search_field[ 'x_to' ]
            y_to = search_field[ 'y_to' ]
            upper_tangent = search_field[ 'upper_tangent' ]
            lower_tangent = search_field[ 'lower_tangent' ]

        # set increments
        if increments:
            n_x = increments[ 'n_x' ]
            n_y = increments[ 'n_y' ]
            n_r = increments[ 'n_r' ]
        
        n_x = max( n_x, 2 )
        n_y = max( n_y, 2 )
        n_r = max( n_r, 2 )

        x_incr = (x_to-x_from) / (n_x-1)
        y_incr = (y_to-y_from) / (n_y-1)
        t_incr = (upper_tangent-lower_tangent) / (n_r-1)

        calc_list = []

        for i in range( n_x ): # prints progression
            x_center = x_from + i*x_incr
            for j in range( n_y ):
                y_center = y_from + j*y_incr
                for k in range( n_r ):
                    radius = y_center - (lower_tangent + k*t_incr)
                    
                    calc_list.append( [ x_center, y_center, radius ] )

#        keys = []
#        for c in calc_list:
#            keys.append( self.fs_manager.gen_fs_key( c[0], c[1], c[2] ) )
#        self.fs_manager.def_fs_dict( set(keys) )

        for i in tqdm( range( len(calc_list) ) ): # prints progression
            c=calc_list[i]
            self.calc_single_circle( c[0], c[1], c[2], clear_FS=False )

        grid = { 
            'x': [ x_from, x_to, x_to, x_from, x_from ], 
            'y': [ y_from, y_from, y_to, y_to, y_from ]
        }

        self.fs_manager.sort( grid )

    def calc_single_circle( self, x_center, y_center, radius, clear_FS=True ):
        if clear_FS:
            self.fs_manager.clear_fs()
        
        self.calc_circle( x_center, y_center, radius )
        fs_id = self.fs_manager.gen_fs_key(x_center, y_center, radius)
        self.populate_lamellas( fs_id )
        self.calc_factor_of_safety( fs_id )


    def calc_factor_of_safety( self, fs_id ):
        fs = self.fs_manager.fs[ fs_id ]
        self.calc_fs( fs )


    def calc_fs( self, fs ):
        eps = 0.0001
        max_it = 100
        k = 0
        F_k = 1 # iterate this one
        F_temp = 0


        while abs(F_k-F_temp) > eps:
            F_temp = F_k # last iteration
            sum_Q_e = 0 # no implemented
            sum_Wi_x_ti = 0 # also called p*delta_x*x_t
            sum_T = 0 # shear force

            lamellas = fs['lamellas'].lamellas
            for i in range( len(lamellas) ):
                lamelle = lamellas[i]
                sum_T += lamelle.get_shear_force( F_k )
                sum_Q_e += 0
                sum_Wi_x_ti += lamelle.weight * (lamelle.x_m_cen-fs['center_x'])

            F_k = ( fs['radius'] * sum_T ) / ( sum_Q_e + sum_Wi_x_ti )

            k += 1
            if k > max_it:
                break
        
        fs['safety'] = F_k

        if False:
            print('sum_T: ' + str(round(sum_T,3)))
            print('sum_Qe: ' + str(round(sum_Q_e,3)))
            print('sum_wtxi: ' + str(round(sum_Wi_x_ti,3)))


    def populate_lamellas( self, fs_id ):
        fs = self.fs_manager.fs[ fs_id ]
        lamellas, undrained = self.calc_lamellas( fs['x'], fs['y'] )
        fs['lamellas'] = lamellas
        fs['undrained'] = undrained


    def calc_lamellas( self, x_list, y_list ):
        lamellas = lam( model_ref=self.model )
        undrained = lamellas.populate( x_list, y_list )
        return [lamellas, undrained]


    def calc_circle( self, x_center, y_center, radius ):
        # get all layer intersections with circle
        layer_icts = self.calc_icts( x_center, y_center, radius ) # calc all circle-layer intersects

        # separate icts by failure surfaces (FS)
        raw_failure_surfaces = self.filter_icts( layer_icts, x_center, y_center, radius ) # find circle intervals penetrating soil, delete icts outside intervals

        # get all layer breakpoints within FS
        raw_fs_x_vals = self.add_intermediate_points( raw_failure_surfaces )

        # calc additional lamelle points
        fs_x_vals = self.add_add_extra_points( raw_fs_x_vals, raw_failure_surfaces, x_center, y_center, radius )
        fs_y_vals = self.calc_fs_y_vals( fs_x_vals, x_center, y_center, radius )

        # save all FS
        self.save_fs( fs_x_vals, fs_y_vals, x_center, y_center, radius )


    def save_fs( self, x_vals, y_vals, x_cen, y_cen, r ):
        for i in range( len( x_vals) ):
            self.fs_manager.add_fs( x_cen, y_cen, r, x_vals[i][0], [x_vals[i], y_vals[i]] )


    def calc_fs_y_vals( self, fs_x_vals, x_center, y_center, radius ):
        y = []
        for x_list in fs_x_vals:
            tmp_y = []
            for x in x_list:
                tmp_y.append( self.calc_circle_y_from_x( x, x_center, y_center, radius ) )
            y.append( tmp_y )    
        return y


    def calc_circle_y_from_x( self, x, x_cen, y_cen, r ):
        return y_cen - (r**2-(x-x_cen)**2)**0.5


    def add_add_extra_points( self, fs_x_vals, raw_FS, x_center, y_center, radius ):
        for i in range( len(raw_FS) ):
            n_points = len( fs_x_vals[i] )
            if n_points < self.n_lamelle: # add when needed, never prune
                n_new = (self.n_lamelle-n_points)+1 # n_points to add

                # calculate FS opening
                x_1, y_1 = [raw_FS[i][0][1][0], raw_FS[i][0][2][0]] # startpoint FS
                x_2, y_2 = [raw_FS[i][0][1][1], raw_FS[i][0][2][1]] # endpoint

                v_1 = [ x_1-x_center, y_1-y_center ]
                v_2 = [ x_2-x_center, y_2-y_center ]

                theta_fs = self.calc_vector_angle( v_1, v_2 )
                start_theta = arctan2( v_1[1], v_1[0] )

                delta_theta = theta_fs / (n_new + 1)

                new_x_vals = self.calc_new_x_vals( radius, x_center, start_theta, delta_theta, n_new )

                combined_x = new_x_vals + fs_x_vals[i]
                combined_x.sort()

                fs_x_vals[i] = combined_x

        return fs_x_vals


    def calc_new_x_vals( self, radius, x_center, start_theta, delta_theta, n_new ):
        new_x_vals = []
        for i in range(n_new):
            new_x_vals.append( radius*cos(delta_theta*(i+1) + start_theta) + x_center )

        return new_x_vals


    def calc_vector_angle( self, v1, v2 ):
        cos_theta = ( (v1[0]*v2[0]) + (v1[1]*v2[1]) ) / ( (v1[0]**2+v1[1]**2)**0.5 * (v2[0]**2+v2[1]**2)**0.5 )
        return arccos(cos_theta)


    def add_intermediate_points( self, raw_fs_es ):
        tmp_fs = copy.deepcopy(raw_fs_es)
        for i in range( len(tmp_fs) ): # failure surface
            for j in range(len(tmp_fs[i])): # layer
                layer_x_vals = self.model.layers[ tmp_fs[i][j][0] ].x

                tmp_x_list = []

                for half_index in range( len(tmp_fs[i][j][1])//2 ): # each dip into layer
                    x_from = tmp_fs[i][j][1][2*half_index]
                    x_to = tmp_fs[i][j][1][2*half_index+1]
                    idx = self.indexes_between_xcoords( layer_x_vals, x_from, x_to, False )
                    
                    tmp_x_list.append( x_from )
                    if idx: # add intermediate points
                        tmp_x_list += [layer_x_vals[k] for k in idx]
                    tmp_x_list.append( x_to )

            tmp_fs[i][j] = tmp_x_list # return x_vals
            
        return tmp_fs[0] ############ check for multiple FS


    def filter_icts( self, layer_icts, x_center, y_center, radius ):
        
        terrain = self.model.layers[0]
        terr_ict_x, terr_ict_y = layer_icts[0][1]

        dipping_icts = [] # index for start of failure surface for a given circle

        if len( terr_ict_x ) > 1:
            for i in range( len( terr_ict_x ) -1 ): # count dips into terrain
                if self.circle_penetrating_terrain( terrain, terr_ict_x[i], terr_ict_y[i], terr_ict_x[i+1], terr_ict_y[i+1], x_center, y_center, radius ):
                    dipping_icts.append( i )

        raw_failure_surfaces = []
        for j in range( len(dipping_icts) ):
            tmp_fs = []
            x_from = layer_icts[0][1][0][dipping_icts[j]] # terrain, coordinates, x, with index j
            x_to = layer_icts[0][1][0][dipping_icts[j]+1] # ... +1

            for l_ict in layer_icts:
                idx = self.indexes_between_xcoords( l_ict[1][0],x_from,x_to, True )
                tmp_fs.append([ l_ict[0], [l_ict[1][0][k] for k in idx], [l_ict[1][1][k] for k in idx] ])
            raw_failure_surfaces.append(tmp_fs)

        return raw_failure_surfaces


    def circle_penetrating_terrain( self, terrain, x_1, y_1, x_2, y_2, x_center, y_center, radius ):
        # find "surface point" on terrain between ict points (from list or coordinate midpoint)
        indexes = self.indexes_between_xcoords( terrain.x, x_1, x_2 )
        if indexes:
            xs_1 = terrain.x[indexes[0]]
            ys_1 = terrain.y[indexes[0]]
        else:
            xs_1 = (x_1+x_2)/2
            ys_1 = (xs_1 - x_1) * (y_2-y_1)/(x_2-x_1) + y_1

        # calc y from circle
        y_circ = y_center - (radius**2-(xs_1-x_center)**2)**0.5

        return ys_1 > y_circ # circle below surface?


    def indexes_between_xcoords( self, x_list, x_1, x_2, include=False ):
        i_list = []
        for i in range(len(x_list)):
            if include:
                if x_list[i]>=x_1 and x_list[i]<=x_2:
                    i_list.append(i)
            else:
                if x_list[i]>x_1 and x_list[i]<x_2:
                    i_list.append(i)
        return i_list


    def calc_icts( self, x_center, y_center, radius ):
        lay_icts = []
        i = 0

        for layer in self.model.layers:
            intersects = self.all_circle_intersects( layer, x_center, y_center, radius )
            lay_icts.append( [i, intersects] )
            i += 1

        return lay_icts


    def all_circle_intersects( self, terrain, x_cen, y_cen, r ):
        xs = []
        ys = []
        for i in range( len( terrain.x ) - 1 ):
            intersects, xsi, ysi = self.calc_intersects( terrain.x[i], terrain.y[i], terrain.x[i+1], terrain.y[i+1], x_cen, y_cen, r )
            if intersects:
                for x, y in zip( xsi, ysi ):
                    xs.append( x )
                    ys.append( y )

        return [ xs, ys ]


    def calc_intersects(self, x_1, y_1, x_2, y_2, x_cen, y_cen, r):
        """Calculates intersect between 
               - a line defined by two points (x_1,y_1) and (x_2, y_2) and 
               - a circle defined by a center point (x_cen,y_cen) with a radius r

               formula solved is:
               M * ( x-x_1 ) + y_1 = ( r^2 - (x-x_cen)^2 ) + y_cen
               where M = ( y_2-y_1 ) / ( x_2-x_1 )

               returns intersects within line segment: (x_1,y_1) - (x_2,y_2)
        """

        xsi = [] # init empty temp vars
        ysi = []

        if x_2 == x_1: # vertical segment
            if x_1 >= x_cen-r and x_1 <= x_cen+r: # within radius of circle
                if x_1 == x_cen-r: # double root
                    xsi = [ x_1, x_1 ]
                    ysi = [ y_cen, y_cen ]
                elif x_1 == x_cen+r: # double root
                    xsi = [ x_1, x_1 ]
                    ysi = [ y_cen, y_cen ]
                else:
                    xsi = [ x_1, x_1 ]
                    ysi = [ y_cen - (r**2-(x_1-x_cen)**2)**0.5, y_cen + (r**2-(x_1-x_cen)**2)**0.5 ]
            else:
                pass # outside ( x_cen ± r ): do nothing
        else:
            M = ( y_2-y_1 ) / ( x_2-x_1 ) # this should now be safe
            a = ( 1 + M**2 ) # never 0!
            b = 2 * ( y_1*M - x_1*M*M - y_cen*M - x_cen )
            c = ( (x_1*M)**2 + y_1**2 + y_cen**2 + x_cen**2 - r**2 ) + 2*( M*( x_1*y_cen - x_1*y_1 ) - y_1*y_cen )

            if (b**2-4*a*c) >= 0: # real roots
                if (b**2-4*a*c) == 0: # perfect square - double root
                    rx_1 = -b/(2*a)
                    ry_1 = M * ( rx_1-x_1 ) + y_1

                    xsi = [ rx_1, rx_1 ]
                    ysi = [ ry_1, ry_1 ]

                else: # two real roots
                    rx_1 = ( -b - (b**2-4*a*c)**0.5 ) / ( 2*a )
                    rx_2 = ( -b + (b**2-4*a*c)**0.5 ) / ( 2*a )
                    ry_1 = M * ( rx_1-x_1 ) + y_1
                    ry_2 = M * ( rx_2-x_1 ) + y_1

                    xsi = [ rx_1, rx_2 ]
                    ysi = [ ry_1, ry_2 ]

            else: # complex roots ( no rel solution )
                pass # do nothing

        xr = [] # init empty return vars
        yr = []
        ict = False

        for x, y in zip( xsi, ysi ): # filter out icts within point range
            if x >= x_1 and x <= x_2 and y >= y_1 and y <= y_2:
                xr.append( x )
                yr.append( y )

        if xr:
            ict = True

        return ict, xr, yr

    def circle_bottom_y_from_x( self, x, x_cen, y_cen, radius ):
        return y_cen - (radius**2-(x-x_cen)**2)**0.5


    class failure_surface_manager(): # failure surface manager
        def __init__( self ):
            self.fs = {} # id & fs-es

            self.grid = None
            self.score_id = {}
        
        def def_fs_dict( self, keys ):
             self.fs = dict.fromkeys( keys )
        
        def sort( self, grid ):
            self.grid = grid
            score = []
            score_id = []

            for key in self.fs:
                score_id.append( key )
                score.append( self.fs[key]['safety'] )

            for val, key in sorted( zip(score,score_id) ):
                self.score_id[key] = val


        def add_fs( self, x_cen, y_cen, radius, start_x, xy ):
            fs_key = self.gen_fs_key( x_cen, y_cen, radius)#, start_x )
            fs_definition = {
                "center_x": x_cen,
                "center_y": y_cen,
                "radius": radius,
                "start_x": start_x,
                "x": xy[0],
                "y": xy[1],
                "safety": None,
                "undrained": None
            }
            self.fs[ fs_key ] = fs_definition
        
        def clear_fs( self ):
            self.fs = {}

        def gen_fs_key( self, x_cen, y_cen, radius): #:, start_x ):
                return ( x_cen, y_cen, radius ) #, start_x )

        def get_all_fs( self ):
            return list( self.fs.values() )

        def get_critical_fs( self ):
            pass