import sqlite3
conn = sqlite3.connect("statparas.db")
cur  = conn.cursor()
exeStr = '''create table StatPara(
                ID          integer primary key,
                NAME        text    NOT NULL,
                SOFTWARE    text,
                FUNCTION    text,
                METHOD      text,
                SOLVENT     text,
                SOURCE      text
            )'''
cur.execute(exeStr)

exeStr = '''create table StatParaData(
                ID          integer primary key autoincrement,
                NAME        text    NOT NULL,
                DIST        text,
                MEAN        float,
                STD         float,
                DEGREE      float,
                STATPARA_ID integer NOT NULL
            )'''
cur.execute(exeStr)
exeStr = '''create table NMRData(
                ID          integer primary key autoincrement,
                NAME        text    NOT NULL,
                C13         float,
                H1          float,
                STATPARA_ID integer NOT NULL
            )'''
cur.execute(exeStr)

exeStr = '''insert into StatPara(ID, NAME, SOFTWARE, FUNCTION, METHOD, SOLVENT, SOURCE) values('%s', '%s', '%s', '%s', '%s', '%s', '%s')'''

cur.execute(exeStr % (1, 'DP4', 'Gaussian09', 'B3LYP', '6-31G(d,p)', 'GasPhase', '10.1021/ja105035r') )
cur.execute(exeStr % (2, 'CP3', 'Gaussian09', 'B3LYP', '6-31G', 'GasPhase', '10.1021/jo900408d') )
exeStr0 = '''insert into StatParaData(NAME, DIST, MEAN, STD, DEGREE, STATPARA_ID) values ('%s', '%s', '%s', '%s', '%s', '%s')'''

tmsStr  = '''insert into NMRData(NAME, C13, H1, STATPARA_ID) values ('%s', '%s', '%s', '%s')'''
import pandas
tb = pandas.read_csv("dp4_plus.csv")
items = pandas.DataFrame(tb, columns=['Item',])
for x in range(1,25):
    item = pandas.DataFrame(tb, columns=[str(x),])
    values = [item.xs(i)[0] for i in range(0,29)]
    cur.execute(exeStr % (x+2, 'DP4+', 'Gaussian09', values[0], values[2], values[1], '10.1021/acs.joc.5b02396'))

    cur.execute(tmsStr % ('TMS', values[3], values[4], x+2))

    for j in range(5, 26, 4):
        cur.execute(exeStr0 % (items.xs(j)[0], 't', values[j], values[j+1], values[j+2], x+2))
        cur.execute(exeStr0 % (items.xs(j)[0], 'n', values[j], values[j+3], 0, x+2))



cur.execute(exeStr0 % ('C13', 't', 0, 1.306, 11.38, 1))
cur.execute(exeStr0 % ('H1', 't', 0, 0.185, 14.18, 1))

cur.execute(exeStr0 % ('correctC13', 'n', 0.547, 0.253, 0, 2))
cur.execute(exeStr0 % ('incorrectC13', 'n', -0.487, 0.533, 0, 2))
cur.execute(exeStr0 % ('correctH1', 'n', 0.478, 0.305, 0, 2))
cur.execute(exeStr0 % ('incorrectH1', 'n', -0.786, 0.835, 0, 2))
cur.execute(exeStr0 % ('correctC13H1', 'n', 0.512, 0.209, 0, 2))
cur.execute(exeStr0 % ('incorrectC13H1', 'n', -0.637, 0.499, 0, 2))

conn.commit()


