from scipy import stats
from CP3 import scaleNMR
import numpy as np

def calculateCDP4( errors, expect, stdev ):
    if not errors: return 1.0
    zs = map( lambda e: abs( ( e - expect ) / stdev ), errors )
    cdp4 = reduce( lambda x, y: x*y, map( lambda z0: stats.norm.cdf( -z0 ), zs ) )
    return cdp4

# t distribution
def calculateTDP4( errors, expect, stdev, degree ):
    if not errors: return 1.0
    zs = map( lambda e: abs( ( e - expect ) / stdev ), errors )
    cdp4 = reduce( lambda x, y: x*y, map( lambda z0: stats.t.cdf( -z0, degree ), zs ) )
    return cdp4

def calculatePDP4( errors, expect, stdev ):
    if not errors: return 1.0
    pdp4 = reduce( lambda x, y: x*y, map( lambda e: stats.norm(expect, stdev).pdf(e), errors ) )
    return pdp4

if __name__ == "__main__": 

    #meanC = 0.0
    #meanH = 0.0
    #stdevC = 2.269372270818724
    #stdevH = 0.18731058105269952
    mean_s = 0.0
    stdev_s = 1.58055
    degree_s = 5.9698
    mean_u_sp = -6.15715
    stdev_u_sp = 2.4881
    degree_u_sp = 6.53169
    mean_u_sp3 = 1.30169
    stdev_u_sp3 = 1.64843
    degree_u_sp3 = 6.28871

    #calc_a_c13 = np.array( [40.56, 27.05, 22.15, 34.93, 58.27, 41.39, 36.06, 41.01, 56.48, 40.18, 33.74, 39.09, 63.03, 39.54, 24.46, 24.25, 41.07 ] )
    #calc_b_c13 = np.array( [40.97, 28.56, 22.53, 37.92, 58.34, 39.27, 34.78, 41.11, 57.31, 41.75, 32.35, 39.35, 63.21, 39.5, 24.43, 24.33, 41.39 ] )
    sp = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0]
    #sp = [0]*len(calc_b_c13)
    calc_a_c13 = np.array([74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000])
    calc_b_c13 = np.array([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
    isomers = [ calc_a_c13, calc_b_c13 ]
    #expValue = np.array( [ 41.0, 26.3, 20.9, 34.6, 56.1, 40.0, 34.5, 41.9, 58.5, 37.4, 32.5, 36.9, 65.1, 40.0, 22.0, 23.0, 43.4] )
    expValue = np.array([76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824])
    #expValue = np.array([73.7292 ,46.5231 ,175.5017, 10.9476 , 141.5293, 125.9367, 128.0500, 127.3075,  60.524, 13.9214])
    cdp4_s = []
    cdp4_u = []
    for isomer in isomers:
        scaledValue = scaleNMR( isomer, expValue )
        scaledErrors = map( lambda x, y: x - y, scaledValue, expValue )
        #cdp4_s.append( calculateCDP4( scaledErrors, mean_s, stdev_s ) )
        cdp4_s.append( calculateTDP4( scaledErrors, mean_s, stdev_s, degree_s ) )


        Errors = map( lambda x, y: x - y, isomer, expValue )

        exp_sp = filter( None, map( lambda x, y: y if x == 1 else None, sp, expValue ))
        errors_sp = filter( None, map( lambda x, y: y if x == 1 else None, sp, Errors ))
        isomer_sp = filter( None, map( lambda x, y: y if x == 1 else None, sp, isomer ))

        errors_sp3 = filter( None, map( lambda x, y: y if x == 0 else None, sp, Errors ))
        exp_sp3 = filter( None, map( lambda x, y: y if x == 0 else None, sp, expValue ))
        isomer_sp3 = filter( None, map( lambda x, y: y if x == 0 else None, sp, isomer ))

        cdp4_u.append( calculateTDP4( errors_sp, mean_u_sp, stdev_u_sp, degree_u_sp )*calculateTDP4( errors_sp3, mean_u_sp3, stdev_u_sp3, degree_u_sp3 ) )

    print "scaled: ", cdp4_s
    print "unscaled: ", cdp4_u
    print map( lambda x: 100*x/sum(cdp4_s), cdp4_s )
    print map( lambda x: 100*x/sum(cdp4_u), cdp4_u )
    cdp4 = map( lambda x, y: x*y, cdp4_s, cdp4_u )
    print map( lambda x: 100*x/sum(cdp4), cdp4 )

