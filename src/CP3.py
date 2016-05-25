import numpy as np
from scipy.stats import linregress

# calculate MAE
def calculateMae( calc, exp ):
    return np.mean( np.abs( calc - exp ) )

# scale the NMR shift
def scaleNMR( calc, exp ):
    slope, intercept, r_value, p_value, std_err = linregress( exp, calc )
    return map( lambda x: ( x - intercept ) / slope, calc )

# calculate MAEprime
def calculateMaePrime( calc, exp ):
    scaled = scaleNMR( calc, exp )
    return np.mean( np.abs( scaled - exp ) )

# calculate cpd
def calculateCpd( calc, exp, calcAv, expAv ):
    deltaExp = exp - expAv
    deltaCalc = calc - calcAv
    return np.sum( deltaExp * deltaCalc ) / np.sum( deltaExp * deltaExp )
    
# calculate cpe
def calculateCpe( calc, exp, calcAv, expAv ):
    deltaExp = exp - expAv
    deltaCalc = calc - calcAv

    denom = np.sum( deltaExp * deltaExp )
    num = np.sum( map( lambda x, y: y**3 / x if abs(x/y) > 1 else y * x, deltaCalc, deltaExp ) )

    return num / denom
# calculate cph
def calculateCph( calc, exp, calcAv, expAv ):
    deltaExp = exp - expAv
    deltaCalc = calc - calcAv

    denom = np.sum( deltaExp * deltaExp )
    num = np.sum( map( lambda x, y: y**3 / x if x/y > 1 else y * x, deltaCalc, deltaExp ) )

    return num / denom

# read data
def readNMRData( filename ):
    data = {}
    with open( filename, 'r' ) as f:
        while 1:
            line = f.readline().strip()
            if line.startswith('C'):
                label, shift = line.split()
                data[label] = float(shift)
            if line.startswith('H'):
                label, shift = line.split()
                data[label] = float(shift)
            #if line.startswith('EXCH'):
            if not line: break
    return data

if __name__  == "__main__":
    exp_b_c13 = np.array([76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824])
    exp_a_c13 = np.array([73.7292 ,46.5231 ,175.5017, 10.9476 , 141.5293, 125.9367, 128.0500, 127.3075,  60.524, 13.9214])
    calc_a_c13 = np.array([74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000])
    calc_b_c13 = np.array([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
    exp_b_1h   = np.array([7.30500, 7.30500, 7.30500, 4.74540, 4.17875, 4.17875, 2.79690, 1.25200, 1.01625])
    calc_a_1h  = np.array([5.27250801496412000000, 2.46511592440598000000, 1.13712631214753000000, 7.43468127418946000000, 7.37704065013850000000, 7.22706900674204000000, 4.15874856882326000000, 1.18811375598228000000, 3.14715999970654000000, 4.16130813951634000000])
    calc_b_1h  = np.array([4.75707855963232000000, 2.74337567543689000000, 1.37752782731981000000, 7.34502691786925000000, 7.34329546303314000000, 7.19890290134792000000, 3.91711685472531000000, 0.90984100630663400000, 4.00544909417348000000, 3.89399303536644000000])

    print "CP1: ", calculateCpd( calc_a_c13, exp_a_c13, calc_b_c13, exp_b_c13 )
    print "CP2: ", calculateCpe( calc_a_c13, exp_a_c13, calc_b_c13, exp_b_c13 )
    print "CP3: ", calculateCph( calc_a_c13, exp_a_c13, calc_b_c13, exp_b_c13 )

    def exchange( array, x, y ):
        array[x], array[y] = array[y], array[x]

    exchange( exp_b_c13, 3, 9 )
    exchange( exp_a_c13, 3, 9 )
    print "CP1: ", calculateCpd( calc_a_c13, exp_b_c13, calc_b_c13, exp_a_c13 )
    print "CP2: ", calculateCpe( calc_a_c13, exp_b_c13, calc_b_c13, exp_a_c13 )
    print "CP3: ", calculateCph( calc_a_c13, exp_b_c13, calc_b_c13, exp_a_c13 )
