from itertools import chain, product, permutations, combinations
from functools import partial
from nmrdata import NMRData, Molecular
from comparemethods import calculateTDP4, calculateCC, calculateMae, calculateCph, getBayersProbabilities
from predata import scaledValue
class StatisticalNMR:
    def __init__(self):
        self.nmrdata = []
        self.molecular = Molecular()

    def readDataFromFile(self, filename):
        raw_calc_data = {}
        raw_exp_data = {}
        with open(filename, "r") as f:
            while True:
                line = f.readline()
                if not line: break
                if line.startswith("#"): continue
                if line.startswith("NAME"):
                    self.molecular.setName( line.lstrip("NAME").strip() )
                    continue
                if line.startswith("DIAS"):
                    items = line.strip().split()
                    self.molecular.setNameofDias(items[1], items[2])
                    continue
                if line.startswith("INDX"):
                    Labels = line.strip().split()[1:]
                    calcItems = filter( lambda line: line.islower(), Labels )
                    expItems  = filter( lambda line: line.isupper(), Labels )
                    continue
                if line.startswith("%BEGIN"):
                    _exchange_list = []
                    _label = []
                    _hybs  = []
                    for item in calcItems:
                        raw_calc_data[item] = []
                    for item in expItems:
                        raw_exp_data[item] = []
                    if not Labels: raise SyntaxError("the format of input file has error!")
                    if( len(line.split()) != 2 ): raise SyntaxError("the format of input file has error!")
                    dtype = line.split()[1]
                    while True:
                        line = f.readline()
                        if line.startswith("%END"): break
                        if line.startswith("EXCH"):
                            items = line.strip().split()[1:]
                            _exchange_list.append(items)
                            continue
                        items = line.strip().split()
                        _label.append(items[0])
                        for i, item in enumerate(Labels):
                            if item.isupper(): raw_exp_data[item].append(items[i+1])
                            elif item.islower(): raw_calc_data[item].append(items[i+1])
                            else: _hybs.append(items[i+1])
                    _nmrdata = NMRData(dtype)
                    _nmrdata.setLabel( _label )
                    self.molecular.setLabel(dtype, _label)
                    self.molecular.setHybs(dtype, _hybs)
                    avg = lambda x: 1.0*sum(x)/len(x)
                    str2float = lambda x: avg(map(float, x.split("-"))) if "-" in x else float(x)
                    for key in raw_calc_data:
                        _nmrdata.setCalcData(key, map(float, raw_calc_data[key]))
                    for key in raw_exp_data:
                        _nmrdata.setExpData(key, map(str2float, raw_exp_data[key]))
                    _nmrdata.setExchangeLabel(_exchange_list)
                    self.nmrdata.append(_nmrdata)
                    continue

    def printNMR(self, label, calc, exp, dtype):
        calc_items = sorted(calc.keys())
        exp_items = sorted(exp.keys())
        data_list = calc_items + exp_items

        label_str = "\t".join([format(item, "<4") for item in ['',] + data_list])
        print("\t"+label_str)

        for i in range(len(label)):
            items = [label[i],] + [format(calc[item][i], ".2f") for item in calc_items ] + \
                [format(exp[item][i], ".2f") for item in exp_items]
            label_str = "\t".join([format(item, "<4") for item in items])
            print("\t"+label_str)

    def checkExchangeItem(self, raw_item ,item):
        exchange_items = []
        for x, y in zip(raw_item, item):
            if x == y: continue
            else: exchange_items.append("%s => %s" % ( x, y ))
        if not exchange_items: return "Raw data"
        return ", ".join( exchange_items )

    def calculateCp3(self):
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

        for _nmrdata in self.nmrdata:
            if(len(_nmrdata.exp_data) < 2):
                raise(Exception("CP3 only for pair dias"))
            exp_per = combinations(_nmrdata.exp_data, 2)
            calc_com = combinations(_nmrdata.calc_data, 2)

            temp = map(lambda x: zip(*x), product(exp_per, calc_com))

            cp3_s = []
            for item in temp:
                argv = []
                item_str = "".join(chain.from_iterable(item))
                pre, av = item
                argv = [_nmrdata.calc_data[pre[1]], _nmrdata.exp_data[pre[0]], 
                        _nmrdata.calc_data[av[1]], _nmrdata.exp_data[av[0]]]
                temp = {}
                temp[item_str] = calculateCph(*argv)

                argv = [_nmrdata.calc_data[pre[1]], _nmrdata.exp_data[av[0]], 
                        _nmrdata.calc_data[av[1]], _nmrdata.exp_data[pre[0]]]
                item_str = "".join((av[0], pre[1], pre[0], av[1]))
                temp[item_str] = calculateCph(*argv)
                cp3_s.append(temp)

            bayers13c = partial( getBayersProbabilities, uv_correct=(expect_c, stdev_c), uv_incorrect=(iexpect_c, istdev_c) )
            bayers1h = partial( getBayersProbabilities, uv_correct=(expect_h, stdev_h), uv_incorrect=(iexpect_h, istdev_h) )
            cp3_p = {}
            if _nmrdata.dtype == "13C":
                for item in cp3_s:
                    key = "/".join(item.keys())
                    cp3_p[key] = bayers13c(*item.values())
                    key = "/".join(item.keys()[::-1])
                    cp3_p[key] = bayers13c(*item.values()[::-1])
            elif _nmrdata.dtype == "1H":
                for item in cp3_s:
                    key = "/".join(item.keys())
                    cp3_p[key] = bayers1h(*item.values())
                    key = "/".join(item.keys()[::-1])
                    cp3_p[key] = bayers1h(*item.values()[::-1])
            print("CP3: " + "; ".join(["{0}:{1:.2f}%".format(x,100*cp3_p[x]) for x in cp3_p]))

    def report(self):
        for _nmrdata in self.nmrdata:
            for i, item in enumerate(_nmrdata.all_exchange_list):
                temp_list = lambda y: [_nmrdata.label[x] for x in chain.from_iterable(y)]
                print("=="*10+"("+str(i+1)+")"+"=="*10)
                print(self.checkExchangeItem(temp_list(_nmrdata.all_exchange_list[0]), temp_list(item)))
                print("--"*20)
                temp_exp_data = _nmrdata.exchangeNMR(item)
                self.printNMR(_nmrdata.label, _nmrdata.calc_data, temp_exp_data, _nmrdata.dtype)
                meanC = 0.0
                meanH = 0.0
                stdevC = 2.306
                stdevH = 0.187
                degreeC = 11.38
                degreeH = 14.18
                for iexp in temp_exp_data:
                    cdp4_s = {}
                    cc_s = {}
                    mae_s = {}
                    for icalc in _nmrdata.calc_data:
                        scaled_value = scaledValue( _nmrdata.calc_data[icalc], temp_exp_data[iexp] )
                        if _nmrdata.dtype == "13C": cdp4_s[icalc] = calculateTDP4(scaled_value, temp_exp_data[iexp], meanC, stdevC, degreeC)
                        if _nmrdata.dtype == "1H":  cdp4_s[icalc] = calculateTDP4(scaled_value, temp_exp_data[iexp], meanH, stdevH, degreeH)
                        cc_s[icalc] = calculateCC(scaled_value, temp_exp_data[iexp])
                        mae_s[icalc] = calculateMae(scaled_value, temp_exp_data[iexp])

                    cdp4_s = map(lambda x: (iexp+x, "{0:.2f}%".format(100*cdp4_s[x]/sum(cdp4_s.values()))),  cdp4_s)
                    cc_s = map(lambda x: (iexp+x, "{0:.6f}".format(cc_s[x])),  cc_s)
                    mae_s = map(lambda x: (iexp+x, "{0:.6f}".format(mae_s[x])),  mae_s)
                    print("DP4: " + "; ".join(["{0}:{1}".format(x,y) for x, y in cdp4_s]))
                    print("CC: " + "; ".join(["{0}:{1}".format(x,y) for x, y in cc_s]))
                    print("MAE: " + "; ".join(["{0}:{1}".format(x,y) for x, y in mae_s]))
                self.calculateCp3()

if __name__ == "__main__":
    import sys
    s = StatisticalNMR()
    filename = sys.argv[1]
    s.readDataFromFile(filename)
    s.report()
    #s.calculateCp3()
