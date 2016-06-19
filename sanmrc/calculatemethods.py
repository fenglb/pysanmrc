from scipy import stats
from functools import partial, reduce
from itertools import compress
from predata import scaledValue
import sqlite3
class StatPara:
    def __init__(self, db_filename):
        self.db_filename = db_filename
        self.conn = sqlite3.connect(self.db_filename)

    def getStatParaID(self, name, soft, function, method, solvent):
        statparaselect_str = '''select ID from StatPara 
                                where NAME=? and 
                                SOFTWARE=? and 
                                FUNCTION=? and 
                                METHOD=? and
                                SOLVENT=?
                                '''
        cursor = self.conn.execute(statparaselect_str, (name, soft, function, method, solvent))
        try:
            return next(cursor)[0]
        except StopIteration:
            return None
    def getStatParaData(self, name, soft, function, method, solvent, dtype, dist):
        statparadataselect_str = '''select MEAN, STD, DEGREE from StatParaData
                                where NAME=? and
                                DIST=? and
                                STATPARA_ID=?
                                '''
        statpara_id = self.getStatParaID(name, soft, function, method, solvent)
        if( not statpara_id ): return None
        cursor = self.conn.execute(statparadataselect_str, (dtype, dist, statpara_id,) )
        try:
            return next(cursor)
        except StopIteration:
            return None
    def getTMSData(self, name, soft, function, method, solvent):
        statpara_id = self.getStatParaID(name, soft, function, method, solvent)
        if( not statpara_id ): return None
        tmsselect_str = '''select C13, H1 from NMRData
                                where NAME = 'TMS' and
                                STATPARA_ID=?
                                '''
        cursor = self.conn.execute(tmsselect_str, (statpara_id,) )
        try:
            return next(cursor)
        except StopIteration:
            return None

# calculate cph
def calculateCph(calc, exp, calcAv, expAv):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise Exception("the length of the exp list must be the same as the calc")
    delta_exp   = (x - y for x, y in zip(exp, expAv))
    delta_calc  = (x - y for x, y in zip(calc, calcAv))

    selector = lambda x, y: y**3 / x if x/y > 1 else y*x

    temp = ((y**2, selector(x,y)) for x, y in zip(delta_calc, delta_exp))
    denom, num = reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), temp )
    return num / denom
# bayers probabilities
def getBayersProbabilities(r1, r2, uv_correct, uv_incorrect):
    expect,  stdev  = uv_correct
    iexpect, istdev = uv_incorrect
    
    probility = lambda x, u, sigma: stats.norm.cdf(-1.0 * abs(x - u) / sigma)

    probility_corr = partial(probility, u=expect, sigma=stdev)
    probility_incorr = partial(probility, u=iexpect, sigma=istdev)

    return probility_corr(r1) * probility_incorr(r2) * 0.5 /  \
          (probility_corr(r1) * probility_incorr(r2) * 0.5 + \
           probility_incorr(r1) * probility_corr(r2) * 0.5 )

# normal distribution
def calculateCDP4(calc, exp, expect, stdev):
    if not ( len(exp) == len(calc) ):
        raise Exception("the length of the exp list must be the same as the calc")

    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.norm.cdf(-1.0 * abs(x - expect) / stdev)

    cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
    return cdp4

# t distribution
def calculateTDP4(calc, exp, expect, stdev, degree):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.t.cdf(-1.0 * abs(x - expect) / stdev, degree)

    cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
    return cdp4

# calculate correlation coefficient
def calculateCC(calc, exp, dtype):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    return stats.pearsonr( exp, calc )[0]
    
# calculate MAE
def calculateMae(calc, exp, dtype):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    return sum( abs(x - y) for x, y in zip(calc, exp) ) / len( calc ) 

def calculateCP3(exp_data_pair, calc_data_pair, uv_correct, uv_incorrect):

    cp3_s = []
    argv = [calc_data_pair[0], exp_data_pair[0], 
            calc_data_pair[1], exp_data_pair[1]]
    cp3_s.append(calculateCph(*argv))

    argv = [calc_data_pair[0], exp_data_pair[1], 
            calc_data_pair[1], exp_data_pair[0]]
    cp3_s.append(calculateCph(*argv))

    bayersp = partial(getBayersProbabilities, uv_correct=uv_correct[:2], uv_incorrect=uv_incorrect[:2] )

    return bayersp(cp3_s[0], cp3_s[1]), bayersp(cp3_s[1], cp3_s[0])

def calculateDP4(calc, exp, dtype, usv):
    mean, stdev, degree = usv
    return calculateTDP4(calc, exp, mean, stdev, degree)

def calculateuTDP4(calc, exp, spX, usv_sp, usv_sp3):

    spX = map(int, spX)
    exp_sp   = [ value for i, value in enumerate(exp)  if spX[i] <= 2 ]
    calc_sp  = [ value for i, value in enumerate(calc) if spX[i] <= 2 ]

    exp_sp3  = [ value for i, value in enumerate(exp)  if spX[i] == 3 ]
    calc_sp3 = [ value for i, value in enumerate(calc) if spX[i] == 3 ]
    
    expect_sp, stdev_sp, degree_sp = usv_sp
    expect_sp3, stdev_sp3, degree_sp3 = usv_sp3

    dp4_sp   = calculateTDP4( calc_sp, exp_sp, expect_sp, stdev_sp, degree_sp )
    dp4_sp3  = calculateTDP4( calc_sp3, exp_sp3, expect_sp3, stdev_sp3, degree_sp3 )
    return dp4_sp * dp4_sp3

def calculateDP4p(calc, exp, spX, usv, usv_sp, usv_sp3):

    uDP4 = calculateuTDP4(calc, exp, spX, usv_sp, usv_sp3 )
    scale_calc = scaledValue(calc, exp)
    sDP4 = calculateDP4(scale_calc, exp, usv)

    return uDP4, uDP4*sDP4

if __name__ == "__main__":
    db_filename = "../data/statparas.db"
    sp = StatPara(db_filename)
    print(sp.getStatParaData('CP3', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'incorrectC13H1', 'n'))
    print(sp.getTMSData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase'))
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'C13', 't'))
    print("="*8)
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'scaledC13', 't'))
    print(sp.getStatParaData('DP4', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'C13', 't'))


    calc_a = [74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000]
    calc_b = ([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
    exp_a = [76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824]
    exp_b = ([73.7292 ,46.5231 ,175.5017, 10.9476 , 141.5293, 125.9367, 128.0500, 127.3075,  60.524, 13.9214])
