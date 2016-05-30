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
