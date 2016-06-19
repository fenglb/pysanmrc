from scipy import stats
from functools import partial, reduce
from itertools import compress
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

class CalculateMethods:
    def __init__(self, db_filename):
        self.db_filename = db_filename
        self.statpara = StatPara(self.db_filename)

        self.DFT_software    = 'Gaussian09'
        self.function        = 'B3LYP'
        self.method          = '6-31G(d,p)'
        self.solvent         = 'GasPhase'

    def setMethod(self, method):
        self.method = method
    def setSolvent(self, solvent):
        self.method = solvent
    def setFunction(self, function):
        self.function = function

    # calculate cph
    def _calculateCph(self, calc, exp, calcAv, expAv):
        if not ( len(exp) == len(calc) == len(calcAv) == len( expAv ) ):
            raise Exception("the length of the exp list must be the same as the calc")
        delta_exp   = (x - y for x, y in zip(exp, expAv))
        delta_calc  = (x - y for x, y in zip(calc, calcAv))

        selector = lambda x, y: y**3 / x if x/y > 1 else y*x

        temp = ((y**2, selector(x,y)) for x, y in zip(delta_calc, delta_exp))
        denom, num = reduce(lambda a, b: (a[0]+b[0], a[1]+b[1]), temp )
        return num / denom
    # bayers probabilities
    def _getBayersProbabilities(self, r1, r2, uv_correct, uv_incorrect):
        expect,  stdev  = uv_correct
        iexpect, istdev = uv_incorrect
        
        probility = lambda x, u, sigma: stats.norm.cdf(-1.0 * abs(x - u) / sigma)

        probility_corr = partial(probility, u=expect, sigma=stdev)
        probility_incorr = partial(probility, u=iexpect, sigma=istdev)

        return probility_corr(r1) * probility_incorr(r2) * 0.5 /  \
              (probility_corr(r1) * probility_incorr(r2) * 0.5 + \
               probility_incorr(r1) * probility_corr(r2) * 0.5 )
    # normal distribution
    def _calculateCDP4(self, calc, exp, expect, stdev):
        if not ( len(exp) == len(calc) ):
            raise Exception("the length of the exp list must be the same as the calc")

        errors = ( x - y for x, y in zip(calc, exp) )
        probility = lambda x: stats.norm.cdf(-1.0 * abs(x - expect) / stdev)

        cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
        return cdp4

    # t distribution
    def _calculateTDP4(self, calc, exp, expect, stdev, degree):
        if len(exp) != len(calc):
            raise Exception("the length of the exp list must be the same as the calc")
        errors = ( x - y for x, y in zip(calc, exp) )
        probility = lambda x: stats.t.cdf(-1.0 * abs(x - expect) / stdev, degree)

        cdp4 = reduce( lambda x, y: x * y, map( probility, errors ) )
        return cdp4
    def _getStatPara(self, name, dtype, dist='t'):
        usv = self.statpara.getStatParaData(name=name, 
                    soft=self.DFT_software, 
                    function=self.function,
                    method=self.method,
                    solvent=self.solvent, 
                    dtype=dtype, dist=dist)
        if not usv:
            raise TypeError("Your Settings donot match for %s. No Statistical paraments found!" % name)
        return usv
    # calculate correlation coefficient
    def calculateCC(self, calc, exp, dtype):
        if len(exp) != len(calc):
            raise Exception("the length of the exp list must be the same as the calc")
        return stats.pearsonr( exp, calc )[0]
        
    # calculate MAE
    def calculateMae(self, calc, exp, dtype):
        if len(exp) != len(calc):
            raise Exception("the length of the exp list must be the same as the calc")
        return sum( abs(x - y) for x, y in zip(calc, exp) ) / len( calc ) 

    def calculateCP3(self, exp_data_pair, calc_data_pair, dtype):

        if dtype=="13C":
            uv_correct = self._getStatPara("CP3", 'correctC13', "n")
            uv_incorrect = self._getStatPara("CP3", 'incorrectC13', "n")
        if dtype=="1H":
            uv_correct = self._getStatPara("CP3", 'correctH1', "n")
            uv_incorrect = self._getStatPara("CP3", 'incorrectH1', "n")
        if dtype=="13C1H":
            uv_correct = self._getStatPara("CP3", 'correctC13H1', "n")
            uv_incorrect = self._getStatPara("CP3", 'incorrectC13H1', "n")

        cp3_s = []
        argv = [calc_data_pair[0], exp_data_pair[0], 
                calc_data_pair[1], exp_data_pair[1]]
        cp3_s.append(self._calculateCph(*argv))

        argv = [calc_data_pair[0], exp_data_pair[1], 
                calc_data_pair[1], exp_data_pair[0]]
        cp3_s.append(self._calculateCph(*argv))

        bayersp = partial(self._getBayersProbabilities, uv_correct=uv_correct[:2], uv_incorrect=uv_incorrect[:2] )

        return bayersp(cp3_s[0], cp3_s[1]), bayersp(cp3_s[1], cp3_s[0])

    def calculateDP4(self, calc, exp, dtype):

        if dtype == "13C":
            mean, stdev, degree = self._getStatPara("DP4", "C13")
        if dtype == "1H":
            mean, stdev, degree = self._getStatPara("DP4", "H1")

        return self._calculateTDP4(calc, exp, mean, stdev, degree)

    def calculateuTDP4(self, calc, exp, spX, usv_sp, usv_sp3):

        spX = map(int, spX)
        exp_sp   = [ value for i, value in enumerate(exp)  if spX[i] <= 2 ]
        calc_sp  = [ value for i, value in enumerate(calc) if spX[i] <= 2 ]

        exp_sp3  = [ value for i, value in enumerate(exp)  if spX[i] == 3 ]
        calc_sp3 = [ value for i, value in enumerate(calc) if spX[i] == 3 ]
        
        expect_sp, stdev_sp, degree_sp = usv_sp
        expect_sp3, stdev_sp3, degree_sp3 = usv_sp3

        dp4_sp   = self._calculateTDP4( calc_sp, exp_sp, expect_sp, stdev_sp, degree_sp )
        dp4_sp3  = self._calculateTDP4( calc_sp3, exp_sp3, expect_sp3, stdev_sp3, degree_sp3 )
        return dp4_sp * dp4_sp3

    def _calculateuDP4(self, calc, exp, dtype, spX):

        if dtype == "13C":
            usv_sp = self._getStatPara("DP4+", "nonscaledC13(sp2)")
            usv_sp3 = self._getStatPara("DP4+", "nonscaledC13(sp3)")
        if dtype == "1H":
            usv_sp = self._getStatPara("DP4+", "nonscaledH1(sp2)")
            usv_sp3 = self._getStatPara("DP4+", "nonscaledH1(sp3)")
        
        return self.calculateuTDP4( calc, exp, spX, usv_sp, usv_sp3 )

    def _calculatesDP4(self, calc, exp, dtype ):

        if dtype == "13C":
            usv = self._getStatPara("DP4+", "scaledC13")
        if dtype == "1H":
            usv = self._getStatPara("DP4+", "scaledH1")
        
        return self._calculateTDP4( calc, exp, *usv )

    def calculateDP4plus(self, calc, exp, dtype, spX):
        
        uDP4 = self._calculateuDP4(calc, exp, dtype, spX)
        sDP4 = self._calculatesDP4(calc, exp, dtype)
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


    calcMethod = CalculateMethods(db_filename)
    calcMethod.setMethod('6-31G(d,p)')
    calc_a = [74.25543650739220000000, 48.44569810240660000000, 176.85339877593200000000, 11.67264309668400000000, 138.21012855822000000000, 121.89729430807700000000, 122.50376880179900000000, 121.05681827872700000000, 58.83075033676600000000, 15.62731956011590000000]
    calc_b = ([76.54606326018590000000, 48.31514100521790000000, 174.56828551337200000000, 18.42221331063060000000, 139.58787236106400000000, 121.50248362594900000000, 122.56072415709800000000, 121.01379975085700000000, 58.34041422828900000000, 15.42125704559360000000])
    exp_a = [76.1218 ,47.0828, 175.6654, 14.2679 , 141.5826, 126.5307, 128.2366, 127.7682, 60.5621, 13.9824]
    exp_b = ([73.7292 ,46.5231 ,175.5017, 10.9476 , 141.5293, 125.9367, 128.0500, 127.3075,  60.524, 13.9214])

    print("Test DP4")
    a = calcMethod.calculateDP4(calc_a, exp_a, '13C')
    b = calcMethod.calculateDP4(calc_b, exp_a, '13C')
    print( 100*a/(a+b), 100*b/(a+b) )
    a = calcMethod.calculateDP4(calc_a, exp_b, '13C')
    b = calcMethod.calculateDP4(calc_b, exp_b, '13C')
    print( 100*a/(a+b), 100*b/(a+b) )

    print("Test sDP4+")
    a = calcMethod._calculatesDP4(calc_a, exp_a, '13C')
    b = calcMethod._calculatesDP4(calc_b, exp_a, '13C')
    print( 100*a/(a+b), 100*b/(a+b) )
    a = calcMethod._calculatesDP4(calc_a, exp_b, '13C')
    b = calcMethod._calculatesDP4(calc_b, exp_b, '13C')
    print( 100*a/(a+b), 100*b/(a+b) )

    print("Test DP4+")
    sp = [3, 3, 2, 3, 2, 2, 2, 2, 3, 3]
    a = calcMethod.calculateDP4plus(calc_a, exp_a, '13C', sp)
    b = calcMethod.calculateDP4plus(calc_b, exp_a, '13C', sp)
    print( "uDP4+: ", 100*a[0]/(a[0]+b[0]), 100*b[0]/(a[0]+b[0]) )
    print( "DP4+: ", 100*a[1]/(a[1]+b[1]), 100*b[1]/(a[1]+b[1]) )
    a = calcMethod.calculateDP4plus(calc_a, exp_b, '13C', sp)
    b = calcMethod.calculateDP4plus(calc_b, exp_b, '13C', sp)
    print( "uDP4+: ", 100*a[0]/(a[0]+b[0]), 100*b[0]/(a[0]+b[0]) )
    print( "DP4+: ", 100*a[1]/(a[1]+b[1]), 100*b[1]/(a[1]+b[1]) )

    print("Test CP3")
    calcMethod.setMethod('6-31G')
    print(calcMethod.calculateCP3([calc_a, calc_b], [exp_a, exp_b], '13C'))
