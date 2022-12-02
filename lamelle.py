"""  
geometry definition

                x2,y2
x1,y1      __--o
      o--‾‾    |
      |        |
      |        |
      |        |
      o--__    |
x4,y4      ‾‾--O
                x3,y3

note vertical boundaries (x1==x4 && x2==x3) 
"""

from numpy import cos, arctan, pi, tan

class lamellas():
    def __init__( self, model_ref):
        self.lamellas = [] # list of classes
        self.model = model_ref


    def populate( self, x_list, y_list ):
        undrained = False
        for i in range( len(x_list)-1 ):
            self.lamellas.append( lamellas.lamelle( self.model, x_list[i], y_list[i], x_list[i+1], y_list[i+1] ) )
            if self.lamellas[ -1 ].undrained:
                undrained = True
        return undrained


    class lamelle():
        def __init__( self, model_ref, x4, y4, x3, y3 ):
            self.soils = []
            
            # define geometry
            self.model = model_ref
            self.x1_4 = x4
            self.y4 = y4
            self.x2_3 = x3
            self.y3 = y3
            self.dx = x3-x4

            self.y1 = self.model.layers[0].y_from_x( x4 ) # top from terrain
            self.y2 = self.model.layers[0].y_from_x( x3 )

            self.x_cen = ( x3 + x4 ) / 2
            self.y_cen = ( y3 + y4 ) / 2

            self.tan_alpha = (y3-y4) / (x3-x4)
            self.alpha = arctan( self.tan_alpha )
            self.cos_alpha = cos( self.alpha )

            self.delta_L = (x3-x4) / self.cos_alpha

            # add soils to lamelle
            self.add_soils()

            # calc porepressure
            self.u = self.calc_pp()

            # calc center of mass
            self.calc_mass_center()

            # get strength params from base soil
            self.get_base_strength()


        def calc_pp( self ):
            if self.model.GW: # calculate PP if GW present
                h_GW = self.model.GW.y_from_x( self.x_cen )
                h_w = ( h_GW-self.y_cen )

                if h_w >= 0:
                    return h_w * self.model.GW.gamma
                return 0
            return 0


        def get_base_strength( self ):
            bottom_soil = self.soils[-1]

            self.gamma = bottom_soil.gamma
            
            self.undrained = bottom_soil.undrained
            self.a = bottom_soil.a
            self.phi = bottom_soil.phi
            self.phi_rad = pi * self.phi / 180
            self.tan_phi = tan( self.phi_rad )
            self.cu = bottom_soil.cu

            self.p = self.weight / self.dx
            self.p_prime = self.p - self.u

            self.T = ( self.p_prime+self.a ) * self.tan_phi * self.dx


        def get_shear_force( self, F ):
            if self.undrained:
                return self.cu * self.delta_L

            # else drained
            self.m_alpha = ( 1+1/F * self.tan_phi * self.tan_alpha ) * self.cos_alpha
            return ( self.T/self.m_alpha )


        def calc_mass_center( self ):
            A = 0
            AW = 0
            x_AW = 0
            y_AW = 0
            for soil in self.soils:
                A += soil.area

                AW += soil.weight
                x_AW += soil.x_m_cen * soil.weight
                y_AW += soil.y_m_cen * soil.weight

            self.area = A
            self.x_m_cen = x_AW/AW
            self.y_m_cen = y_AW/AW
            self.weight = AW


        def get_xy( self ): # for plotting
            x = [ self.x1_4, self.x2_3, self.x2_3, self.x1_4, self.x1_4 ]
            y = [ self.y1, self.y2, self.y3, self.y4, self.y1 ]
            return x, y


        def add_soils( self ):
            bottom_i = self.model.get_layer_index_by_xy( self.x_cen, self.y_cen )

            for i in range( bottom_i+1 ):
                y1 = self.model.layers[i].y_from_x( self.x1_4 )
                y2 = self.model.layers[i].y_from_x( self.x2_3 )

                if i == bottom_i:
                    y3 = self.y3
                    y4 = self.y4
                else:
                    y3 = self.model.layers[i+1].y_from_x( self.x2_3 )
                    y4 = self.model.layers[i+1].y_from_x( self.x1_4 )

                self.soils.append( lamellas.soil_volume( self.model, i, self.x1_4, y1, self.x2_3, y2, self.x2_3, y3, self.x1_4, y4) )


    class soil_volume():
        def __init__( self, model_ref, layer_index, x1, y1, x2, y2, x3, y3, x4, y4 ):            
            self.model = model_ref
            self.layer_index = layer_index
            self.x = [ x1, x2, x3, x4 ]
            self.y = [ y1, y2, y3, y4 ]

            self.area, self.x_m_cen, self.y_m_cen = self.calc_mass_center()
            self.gamma, self.undrained, self.a, self.phi, self.cu = self.get_soil_params( self.layer_index )
            self.weight = self.gamma * self.area


        def get_soil_params( self, layer_index ):
            layer = self.model.layers[ layer_index ]

            return [ layer.gamma, layer.undrained, layer.a, layer.phi, layer.cu ]


        def calc_mass_center( self ):
            x = self.x + [ self.x[0] ] # closed polygon
            y = self.y + [ self.y[0] ]

            # area
            first_sum =0
            for i in range( len(x)-1 ):
                first_sum += ( x[i]*y[i+1] - x[i+1]*y[i] )
            area = -first_sum / 2

            # centroid x
            second_sum = 0
            for i in range( len(x)-1 ):
                second_sum += ( (x[i]+x[i+1]) * (x[i]*y[i+1] - x[i+1]*y[i]) )
            x_m_cen = -1/( 6*area ) * second_sum

            # centroid y
            third_sum = 0
            for i in range( len(x)-1 ):
                third_sum += ( (y[i]+y[i+1]) * (x[i]*y[i+1] - x[i+1]*y[i]) )
            y_m_cen = -1/( 6*area ) * third_sum

            return [area, x_m_cen, y_m_cen]