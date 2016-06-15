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

        getpara = partial(self.statpara.getStatParaData, name='CP3', 
                soft=self.DFT_software, 
                function=self.function,
                method=self.method,
                solvent=self.solvent, 
                dist='n'):
        if dtype=="13C":
            uv_correct = getpara('correctC13')
            uv_incorrect = getpara('incorrectC13')
        if dtype=="1H":
            uv_correct = getpara('correctH1')
            uv_incorrect = getpara('incorrectH1')
        if dtype=="13C1H":
            uv_correct = getpara('correctC13H1')
            uv_incorrect = getpara('incorrectC13H1')
                        )
        cp3_s = []
        argv = [calc_data_pair[0], exp_data_pair[0], 
                calc_data_pair[1], exp_data_pair[1]]
        cp3_s.append(calculateCph(*argv))

        argv = [calc_data_pair[0], exp_data_pair[1], 
                calc_data_pair[1], exp_data_pair[0]]
        cp3_s.append(calculateCph(*argv))

        bayersp = partial(self._getBayersProbabilities, uv_correct[:2], uv_incorrect[:2] )

        return bayersp(cp3_s[0], cp3_s[1]), bayersp(cp3_s[1], cp3_s[0])

    def calculateDP4(self, calc, exp, dtype):
        meanC = 0.0
        meanH = 0.0
        stdevC = 2.306

        stdevH = 0.187
        degreeC = 11.38
        degreeH = 14.18
        if dtype == "13C": return self._calculateTDP4(calc, exp, meanC, stdevC, degreeC)
        if dtype == "1H":  return self._calculateTDP4(calc, exp, meanH, stdevH, degreeH)


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
            usv_sp = [-6.16, 2.49, 6.53]
            usv_sp3 = [1.3, 1.65, 6.29]
        if dtype == "1H":
            usv_sp = [-0.17, 0.19, 13.81]
            usv_sp3 = [-0.06, 0.17, 6.95]
        
        return self._calculateuTDP4( calc, exp, spX, usv_sp, usv_sp3 )

    def _calculatesDP4(self, calc, exp, dtype ):

        if dtype == "13C":
            usv = [-6.16, 2.49, 6.53]
        if dtype == "1H":
            usv = [-0.17, 0.19, 13.81]
        
        return self._calculateTDP4( calc, exp, *usv )

    def calculateDP4plus(self, calc, exp, dtype, spX):
        
        uDP4 = self._calculateuDP4(calc, exp, dtype, spX)
        sDP4 = self._calculatesDP4(calc, exp, dtype)
        return uDP4, uDP4*sDP4

if __name__ == "__main__":
    sp = StatPara("../data/statparas.db")
    print(sp.getStatParaData('CP3', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'incorrectC13H1', 'n'))
    print(sp.getTMSData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase'))
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'C13', 't'))
    print("="*8)
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'scaledC13', 't'))
    print(sp.getStatParaData('DP4', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'C13', 't'))

