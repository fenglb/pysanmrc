from scipy import stats
from functools import partial, reduce
from itertools import compress
# calculate correlation coefficient
def calculateCC( calc, exp ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")
    return stats.pearsonr( exp, calc )[0]
    
# calculate MAE
def calculateMae( calc, exp ):
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

def calculatePDP4( calc, exp, expect, stdev ):
    if len(exp) != len(calc):
        raise Exception("the length of the exp list must be the same as the calc")

    errors = ( x - y for x, y in zip(calc, exp) )
    probility = lambda x: stats.norm(expect, stdev).pdf(x)

    pdp4 = reduce( lambda x, y: probility(x) * probility(y), errors )
    return pdp4

def calculateDP4plus( calc, exp, spX, usv_sp, usv_sp3, distribution="t" ):
    if not ( len(exp) == len(calc) == len(spX) ):
        raise Exception("the length of the exp list must be the same as the calc, and spX")
    if not ( isinstance(usv_sp, list) and isinstance(usv_sp3, list ) ):
        raise Exception("the mean, stdev and degree must be included in a list")
    if not ( len(usv_sp3) >= 2 and len(usv_sp) >= 2 ):
        raise Exception("usv_sp and usv_sp3 must have the mean and stdev if you use normal distribution")
    if distribution == "t" and not ( len(usv_sp3) == 3 and len(usv_sp) == 3 ):
        raise Exception("usv_sp and usv_sp3 must have the mean, stdev and degree if you use t-distribution")

    exp_sp  = [ value for i, value in enumerate(exp) if spX[i] == 1 ]
    calc_sp = [ value for i, value in enumerate(calc) if spX[i] == 1 ]

    exp_sp3  = [ value for i, value in enumerate(exp) if spX[i] == 0 ]
    calc_sp3 = [ value for i, value in enumerate(calc) if spX[i] == 0 ]
    
    expect_sp, stdev_sp, degree_sp = usv_sp
    expect_sp3, stdev_sp3, degree_sp3 = usv_sp3

    dp4_sp  = 0.0
    dp4_sp3 = 0.0

    if( distribution == "t" ):
        dp4_sp   = calculateTDP4( calc_sp, exp_sp, expect_sp, stdev_sp, degree_sp )
        dp4_sp3  = calculateTDP4( calc_sp3, exp_sp3, expect_sp3, stdev_sp3, degree_sp3 )
    elif ( distribution == "c" ):
        dp4_sp   = calculateCDP4( calc_sp, exp_sp, expect_sp, stdev_sp )
        dp4_sp3  = calculateCDP4( calc_sp3, exp_sp3, expect_sp3, stdev_sp3 )
    elif ( distribution == "p" ):
        dp4_sp   = calculatePDP4( calc_sp, exp_sp, expect_sp, stdev_sp )
        dp4_sp3  = calculatePDP4( calc_sp3, exp_sp3, expect_sp3, stdev_sp3 )

    return dp4_sp * dp4_sp3
