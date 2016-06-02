from scipy import stats

def scaledValue( calc, exp ):
    if ( len(calc) != len(exp) ):
        raise("the length of the exp list must be the same as the calc")
    slope, intercept, r_value, p_value, std_err = stats.linregress( exp, calc )
    return list( ( x - intercept ) / slope for x in calc )
