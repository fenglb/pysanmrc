from scipy import stats
# calculate correlation coefficient
def calculateCC( calc, exp ):
    if len(exp) != len(calc):
        raise("the length of the exp list must be the same as the calc")
    return stats.pearsonr( exp, calc )[0]
    
# calculate MAE
def calculateMae( calc, exp ):
    if len(exp) != len(calc):
        raise("the length of the exp list must be the same as the calc")
    return sum( map( lambda x, y: abs(x - y ), calc, exp ) ) / len( calc ) 

# calculate cpd
def calculateCpd( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise("the length of the exp list must be the same as the calc")
    delta_exp   = map( lambda x, y: x - y, exp, expAv )
    delta_calc  = map( lambda x, y: x - y, calc, calcAv )
    return sum( map( lambda x, y: x * y, delta_exp, delta_calc ) ) / sum( map ( lambda x, y: x * y, delta_exp, delta_exp ) )

# calculate cpe
def calculateCpe( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise("the length of the exp list must be the same as the calc")
    delta_exp   = map( lambda x, y: x - y, exp, expAv )
    delta_calc  = map( lambda x, y: x - y, calc, calcAv )
    
    denom       = sum( map( lambda x, y: x * y, delta_exp, delta_exp ) )
    num         = sum( map( lambda x, y: y**3 / x if abs(x/y) > 1 else y * x, delta_calc, delta_exp ) )

    return num / denom

# calculate cph
def calculateCph( calc, exp, calcAv, expAv ):
    if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
        raise("the length of the exp list must be the same as the calc")
    delta_exp   = map( lambda x, y: x - y, exp, expAv )
    delta_calc  = map( lambda x, y: x - y, calc, calcAv )
    
    denom       = sum( map( lambda x, y: x * y, delta_exp, delta_exp ) )
    num         = sum( map( lambda x, y: y**3 / x if x/y > 1 else y * x, delta_calc, delta_exp ) )

    return num / denom

# normal distribution
def calculateCDP4( calc, exp, expect, stdev ):
    if not ( len(exp) == len(calc) ):
        raise("the length of the exp list must be the same as the calc")

    errors = map( lambda x, y: x - y, calc, exp )
    zs = map( lambda e: abs( ( e - expect ) / stdev ), errors )
    cdp4 = reduce( lambda x, y: x*y, map( lambda z0: stats.norm.cdf( -z0 ), zs ) )
    return cdp4

# t distribution
def calculateTDP4( calc, exp, expect, stdev, degree ):
    if not ( len(exp) == len(calc) ):
        raise("the length of the exp list must be the same as the calc")

    errors = map( lambda x, y: x - y, calc, exp )
    zs = map( lambda e: abs( ( e - expect ) / stdev ), errors )
    cdp4 = reduce( lambda x, y: x*y, map( lambda z0: stats.t.cdf( -z0, degree ), zs ) )
    return cdp4

def calculatePDP4( calc, exp, expect, stdev ):
    if not ( len(exp) == len(calc) ):
        raise("the length of the exp list must be the same as the calc")

    errors = map( lambda x, y: x - y, calc, exp )
    pdp4 = reduce( lambda x, y: x*y, map( lambda e: stats.norm(expect, stdev).pdf(e), errors ) )
    return pdp4

def calculateDP4plus( calc, exp, spX, usv_sp, usv_sp3, distribution="t" ):
    if not ( len(exp) == len(calc) == len(spX) ):
        raise("the length of the exp list must be the same as the calc, and spX")
    if not ( isinstance(usv_sp, list) and isinstance(usv_sp3, list ) ):
        raise("the mean, stdev and degree must be included in a list")
    if not ( len(usv_sp3) >= 2 and len(usv_sp) >= 2 ):
        raise("usv_sp and usv_sp3 must have the mean and stdev if you use normal distribution")
    if distribution == "t" and not ( len(usv_sp3) == 3 and len(usv_sp) == 3 ):
        raise("usv_sp and usv_sp3 must have the mean, stdev and degree if you use t-distribution")

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
