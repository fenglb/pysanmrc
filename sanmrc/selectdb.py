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

if __name__ == "__main__":
    sp = StatPara("../data/statparas.db")
    print(sp.getStatParaData('CP3', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'incorrectC13', 'n'))
    print(sp.getTMSData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase'))
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', 'C13', 't'))
    print("="*8)
    print(sp.getStatParaData('DP4+', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'scaledC13', 't'))
    print(sp.getStatParaData('DP4', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', 'C13', 't'))

