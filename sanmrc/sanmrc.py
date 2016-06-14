from itertools import chain, product, permutations, combinations
from functools import partial
from nmrdata import NMRData, Molecular
from comparemethods import calculateDP4, calculateCC, calculateMae, calculateCP3
from predata import scaledValue
class StatisticalNMR:
    def __init__(self):
        self.nmrdata = []
        self.molecular = Molecular()

        map(lambda func: setattr(self, func.__name__, func), 
                [calculateDP4, calculateCC, calculateMae, calculateCP3])
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

    def productPairData(self, data, num):
        for idata in combinations(data, num):
            yield "".join(idata), [data[key] for key in idata]
       
    def productRawData(self):
        for _nmrdata in self.nmrdata:
            for i, item in enumerate(_nmrdata.all_exchange_list):
                temp_list = lambda y: [_nmrdata.label[x] for x in chain.from_iterable(y)]
                print(self.checkExchangeItem(temp_list(_nmrdata.all_exchange_list[0]), temp_list(item)))
                temp_exp_data = _nmrdata.exchangeNMR(item)
                yield _nmrdata.dtype, _nmrdata.label, temp_exp_data, _nmrdata.calc_data

    def calc(self, exps, calcs, funcs, dtype=None):
        results = {}
        for func_name in funcs:
            results[func_name] = {}
            func = getattr(self, func_name, None)
            if not func: continue

            if func_name == 'calculateCP3':
                for exp_label, exp_data in self.productPairData(exps,2):
                    results[func_name][exp_label] = {}
                    for calc_label, calc_data in self.productPairData(calcs,2):
                        scaled_calc_data = map( scaledValue, calc_data, exp_data )
                        results[func_name][exp_label][calc_label] = func(scaled_calc_data, exp_data, dtype)
            else:
                for tab_exp in exps:
                    results[func_name][tab_exp] = {}
                    for tab_calc in calcs:
                        # pre treatment
                        calc = calcs[tab_calc]; exp = exps[tab_exp]
                        scaled_calc_value = scaledValue(calc, exp)

                        results[func_name][tab_exp][tab_calc] = func(scaled_calc_value, exp, dtype)
        return results

    def report(self):
        for dtype, label, exp_data, calc_data in self.productRawData():
            self.printNMR(label, calc_data, exp_data, dtype)

            results = self.calc(exp_data, calc_data, 
                    ['calculateDP4','calculateCP3', 'calculateCC', 'calculateMae'],
                    dtype)
            print results

            cdp4_s = results['calculateDP4']
            #cdp4_s = map(lambda x: (x, "{0:.2f}%".format(100*cdp4_s[x]/sum(cdp4_s.values()))),  cdp4_s)
            #cc_s = map(lambda x: (x, "{0:.6f}".format(cc_s[x])),  cc_s)
            #mae_s = map(lambda x: (x, "{0:.6f}".format(mae_s[x])),  mae_s)
            #print("DP4: " + "; ".join(["{0}:{1}".format(x,y) for x, y in cdp4_s]))
            #print("CC: " + "; ".join(["{0}:{1}".format(x,y) for x, y in cc_s]))
            #print("MAE: " + "; ".join(["{0}:{1}".format(x,y) for x, y in mae_s]))
            #self.calculateCp3()

if __name__ == "__main__":
    import sys
    s = StatisticalNMR()
    filename = sys.argv[1]
    s.readDataFromFile(filename)
    s.report()
    #s.calculateCp3()
