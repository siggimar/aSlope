'''
Possible TODO-s
    # DONE - redo plotting procedure for heatmap (search zone),
        DONE - 1. take known points & use kriging/RBF to fill in unknown areas
    # separate failure surface calculations (geometry) and factor of safety calculations (geotechnics)
    # implement calculation method as separate classes ('OMS', 'Bishop', 'Janbu', ...)
    # implement terrain loads
    # implent varying Su with depth
    # implement ADP
    # fix issues when shear surface hits model boundaries
    # fix issues regarding inclomplete models (extend layers to boundaries)
    # fix geometric issue when circle center is below top of slope

    # automate search for critical fs:
      1. grid calculations (very coarse - overview )
      2. upscale w. meshgrid & kriging
      3. select cutoff FS for detecting lowest FS
      4. grid calculations (medium grid) to rate candidates
      5. grid calculations (fine grid on best candidate)
      6. plot best results, optionally with rest faded in background

    # ...the list goes on
'''

from model import SOIL_MODEL as M
from renderer import MODEL_RENDERER as R
from stability import sircular_cylindric_fs as FS


def main():
    # init soil model, calc_method and renderer
    model = M()
    fs = FS( model, n_lamelle=30 )
    renderer = R( model, fs )


# verification examples (results: [ from Excel ]/[ from this program ]
    # undrained, 2 sircles not critical
    if False: # FØ1 BYGT2001:2022 U43, expected result from Excel with ~x_t in comments /result from here
        model.simple_geom( H=10, L=20, D_ROCK=10, gamma=19, a=10, phi=29, cu=32, undrained=True ) 
        fs.calc_single_circle( 9.03, 20.09, 25.09) # 1a: 1.064/1.063 - OK
        fs.calc_single_circle( 11.57, 10.43, 20.43, clear_FS=False ) # 1b: 1.060/1.059 - OK

    # drained, 2 sircles not critical
    elif False: # FØ2 BYGT2001:2022 U43 - deep circles - GW as in example
        model.simple_geom( H=10, L=20, D_ROCK=10, gamma=19, a=10, phi=29, cu=32, undrained=False )
        model.set_gw( [-12.1576,0,2.5,7.5,12.5,17.5,23,28,30.3490,39.5625],[0,0,0.2474,0.6829,1.0396,1.3177,1.5331,1.6468,1.6733,1.7433] )
        fs.calc_single_circle( 9.03, 20.09, 25.09 ) # 1a: 1.958/1.957 - OK
        fs.calc_single_circle( 11.57, 10.43, 20.43, clear_FS=False ) # 1b: 2.502/2.499 - OK

    # drained, 2 sircles not critical
    elif False: # FØ2 BYGT2001:2022 U43 - deep circles - GW at terrain
        model.simple_geom( H=10, L=20, D_ROCK=10, gamma=19, a=10, phi=29, cu=32, undrained=False )
        model.set_gw( [-15,0,20,40],[0,0,10,10] )
        fs.calc_single_circle( 9.03, 20.09, 25.09 ) # 1a: 1.241/1.228 - OK
        fs.calc_single_circle( 11.57, 10.43, 20.43, clear_FS=False ) # 1b: 1.775/1.773 - OK

    # drained, 2 sircles not critical
    elif False: # FØ2b BYGT2001:2022 U43 - shallow circles
        model.simple_geom( H=10, L=20, D_ROCK=10, gamma=19, a=10, phi=29, cu=32, undrained=False )
        model.set_gw( [-12.1576,0,2.5,7.5,12.5,17.5,23,28,30.3490,39.5625],[0,0,0.2474,0.6829,1.0396,1.3177,1.5331,1.6468,1.6733,1.7433] )
        fs.calc_single_circle( 3.12, 24.27, 24.47 ) # 1a: 1.604/1.604 - OK
        fs.calc_single_circle( 6.64, 21.54, 22.54, clear_FS=False ) # 1b: 1.727/1.726 - OK

    # grid search
    elif True:
        model.simple_geom( H=10, L=20, D_ROCK=10, gamma=19, a=10, phi=29, cu=32, undrained=True )
        #model.set_gw( [50.0,0,2.5,7.5,12.5,17.5,23,28,30.3490,70],[0,0,0.2474,0.6829,1.0396,1.3177,1.5331,1.6468,1.6733,1.7433] ) # fails

        box = { # box, and lowest point bounds
            'x_from': -1,
            'y_from': 11,
            'x_to': 21,
            'y_to': 36,
            'upper_tangent': -0.5,
            'lower_tangent': -10
        }

        inc = {
            'n_x': 20,
            'n_y': 20,
            'n_r': 15
        }

        fs.grid_search( search_field=box, increments=inc )

    if False:# prints lamelle geotemtrical data
        fs.print_geom()
        fs.print_gv()


    renderer.render()

if __name__=='__main__':
    main()