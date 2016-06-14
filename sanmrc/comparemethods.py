from scipy import stats
from functools import partial, reduce
from itertools import compress

# calculate correlation coefficient
def calculateCC( calc, exp, dtype ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    return stats.pearsonr( exp, calc )[0]
    
# calculate MAE
def calculateMae( calc, exp, dtype ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    return sum( abs(x - y) for x, y in zip(calc, exp) ) / len( calc ) 

# calculate cpd
def calculateCpd( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise Exception("the length of the exp list must be the same as the calc")
    delta_exp   = (x - y for x, y in zip(exp, expAv))
    delta_calc  = (x - y for x, y in zip(calc, calcAv))

    selector = lambda x, y: x * y
    temp = ((y**2, selector(x,y)) for x, y in zip(delta_calc, delta_exp))
    denom, num = reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), temp )
    return num / denom

# calculate cpe
def calculateCpe( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise Exception("the length of the exp list must be the same as the calc")
    delta_exp   = (x - y for x, y in zip(exp, expAv))
    delta_calc  = (x - y for x, y in zip(calc, calcAv))
    
    selector = lambda x, y: y**3 / x if abs(x/y) > 1 else y*x
    temp = ((y**2, selector(x,y)) for x, y in zip(delta_calc, delta_exp))
    denom, num = reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), temp )
    return num / denom

# calculate cph
def calculateCph( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise Exception("the length of the exp list must be the same as the calc")
    delta_exp   = (x - y for x, y in zip(exp, expAv))
    delta_calc  = (x - y for x, y in zip(calc, calcAv))

    selector = lambda x, y: y**3 / x if x/y > 1 else y*x

    temp = ((y**2, selector(x,y)) for x, y in zip(delta_calc, delta_exp))
    denom, num = reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), temp )
    return num / denom

# bayers probabilities
def getBayersProbabilities( r1, r2, uv_correct, uv_incorrect ):
    expect,  stdev  = uv_correct
    iexpect, istdev = uv_incorrect
    
    probility = lambda x, u, sigma: stats.norm.cdf(-1.0 * abs(x - u) / sigma)

    probility_corr = partial(probility, u=expect, sigma=stdev)
    probility_incorr = partial(probility, u=iexpect, sigma=istdev)

    return probility_corr(r1) * probility_incorr(r2) * 0.5 /  \
          (probility_corr(r1) * probility_incorr(r2) * 0.5 + \
           probility_incorr(r1) * probility_corr(r2) * 0.5 )

def calculateCP3(exp_data_pair, calc_data_pair, dtype):
    expect_h = 0.478
    stdev_h = 0.303
    iexpect_h = -0.786
    istdev_h = 0.835

    expect_c = 0.547
    stdev_c = 0.253
    iexpect_c = -0.487
    istdev_c = 0.533

    expect_m = 0.512
    stdev_m = 0.209
    iexpect_m = -0.637
    istdev_m = 0.499

    cp3_s = []
    argv = [calc_data_pair[0], exp_data_pair[0], 
            calc_data_pair[1], exp_data_pair[1]]
    cp3_s.append(calculateCph(*argv))

    argv = [calc_data_pair[0], exp_data_pair[1], 
            calc_data_pair[1], exp_data_pair[0]]
    cp3_s.append(calculateCph(*argv))

    bayers13c = partial( getBayersProbabilities, uv_correct=(expect_c, stdev_c), uv_incorrect=(iexpect_c, istdev_c) )
    bayers1h = partial( getBayersProbabilities, uv_correct=(expect_h, stdev_h), uv_incorrect=(iexpect_h, istdev_h) )
    cp3_p = {}
    if dtype == "13C":
        return bayers13c(cp3_s[0], cp3_s[1]), bayers13c(cp3_s[1], cp3_s[0])
    elif dtype == "1H":
        return bayers1h(cp3_s[0], cp3_s[1]), bayers1h(cp3_s[1], cp3_s[0])

# normal distribution
def calculateCDP4( calc, exp, expect, stdev ):
    if not ( len(exp) == len(calc) ):
        raise Exception("the length of the exp list must be the same as the calc")

    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.norm.cdf(-1.0 * abs(x - expect) / stdev)

    cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
    return cdp4

# t distribution
def calculateTDP4( calc, exp, expect, stdev, degree ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.t.cdf(-1.0 * abs(x - expect) / stdev, degree)

    cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
    return cdp4

def calculateDP4(calc, exp, dtype):
    meanC = 0.0
    meanH = 0.0
    stdevC = 2.306
    stdevH = 0.187
    degreeC = 11.38
    degreeH = 14.18
    if dtype == "13C": return calculateTDP4(calc, exp, meanC, stdevC, degreeC)
    if dtype == "1H":  return calculateTDP4(calc, exp, meanH, stdevH, degreeH)


def calculatePDP4( calc, exp, expect, stdev ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")

    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.norm(expect, stdev).pdf(x)

    pdp4 = reduce( lambda x, y: probility(x) * probility(y), errors )
    return pdp4

def calculateuTDP4( calc, exp, spX, usv_sp, usv_sp3 ):

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

def calculateuDP4(calc, exp, dtype, spX):

    if dtype == "13C":
        usv_sp = [-6.16, 2.49, 6.53]
        usv_sp3 = [1.3, 1.65, 6.29]
    if dtype == "1H":
        usv_sp = [-0.17, 0.19, 13.81]
        usv_sp3 = [-0.06, 0.17, 6.95]
    
    return calculateuTDP4( calc, exp, spX, usv_sp, usv_sp3 )
