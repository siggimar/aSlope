class SOIL_MODEL():
    def __init__( self ):
        ''' creates a single layer geometry '''
        self.dy = 0.01 # m
        self.layers = [] # list of (non-intersecting) LAYERS in descending order
        self.GW = None # single LAYER with x, y + gamma_w
        self.rock = None # single LAYER with x, y + high cu
        self.bounds = [ [0, 0], [0, 0] ] # [[x_min, x_max], [y_min, y_max]]

        self.failure_surfaces = [] # list of lists: [ [xs],[ys],Factor of safety ]
        self.sircles = [] # list of circles to plot [ x_cen, y_cen, r ]

        self.layer_color_max_index = 10  # 10: 0-9
        self.layer_color_index = 0


    def simple_geom( self, H, L, D_ROCK, gamma=20, a=0, phi=30, cu=0, undrained=True ):
        ''' creates a simple slope and flat bedrock '''
        n_d = 5
        x = [ -D_ROCK*n_d, 0, L, L+D_ROCK*n_d ]
        y = [ 0, 0, H, H ]
        self.add_layer( x, y, gamma, a, phi, cu, undrained=undrained )
        self.set_rock( x, [ -D_ROCK, -D_ROCK, -D_ROCK, -D_ROCK ] )


    def set_gw( self, x_list, y_list ):
        self.GW = SOIL_MODEL.LAYER( x_list, y_list, gamma=10, cu=0 )


    def add_layer( self, x, y, gamma, a, phi, cu, undrained=True ):
        c_i = self.layer_color_index%self.layer_color_max_index
        self.layer_color_index += 1

        self.layers.append( SOIL_MODEL.LAYER( x, y, gamma, a, phi, cu, undrained=undrained, color_index=c_i ) )
        self.check_bounds( x, y )


    def set_rock( self, x, y ):
        self.rock = SOIL_MODEL.LAYER( x, y, gamma=0.01, cu=999999 )
        self.check_bounds( x, y )


    def set_GW( self, x, y ):
        self.GW = SOIL_MODEL.LAYER( x, y )
        self.check_bounds( x, y )


    def check_bounds( self, x, y ):
        if min(x)<self.bounds[0][0]: self.bounds[0][0]=min(x)
        if max(x)>self.bounds[0][1]: self.bounds[0][1]=max(x)
        if min(y)<self.bounds[1][0]: self.bounds[1][0]=min(y)
        if max(y)>self.bounds[1][1]: self.bounds[1][1]=max(y)


    def get_layer_index_by_xy( self, x_cen, y_cen ):
        for i in range( len(self.layers) ):
            y_cen_layer = self.layers[i].y_from_x( x_cen )
            if  y_cen_layer < y_cen:
                return (i-1)
        return len(self.layers)-1


    def sort_layers( self ):
        # implement later
        pass


    class LAYER():
        def __init__( self, x=[], y=[], gamma=20.0, a=0.0, phi=0.0, cu=0.0, undrained=True, color_index=0 ):
            self.x = x
            self.y = y

            self.undrained = undrained

            self.gamma = gamma

            self.a = a
            self.phi = phi
            self.cu = cu

            self.color_i = color_index
        
        def y_from_x( self, x ):
            i = self.get_low_x_index(x)
            x1 =  self.x[i]
            x2 =  self.x[i+1]
            y1 =  self.y[i]
            y2 =  self.y[i+1]
            
            return (y2-y1) / (x2-x1) * (x-x1) + y1

        def get_low_x_index( self, x ):
            for i in range(len(self.x)-1):
                if self.x[i]<=x and self.x[i+1]>=x:
                    return i
            return (len(self.x)-1)

        def sort_coords( self ): # ensure coords are left to right
            pass