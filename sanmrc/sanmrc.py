from nmrdata import NMRData, Molecular
class StatisticalNMR:
    def __init__(self):
        self.nmrdata = []
        self.molecular = Molecular()

    def readDataFromFile(self, filename, dtype):
        _exchange_list = []
        raw_calc_data = {}
        raw_exp_data = {}
        _label = []
        _hybs  = []
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
                    for item in calcItems:
                        raw_calc_data[item] = []
                    for item in expItems:
                        raw_exp_data[item] = []
                    continue
                if line.startswith("%BEGIN"):
                    if not Labels: raise SyntaxError("the format of input file has error!")
                    while True:
                        line = f.readline()
                        if line.startswith("%END"): break
                        items = line.strip().split()
                        _label.append(items[0])
                        for i, item in enumerate(Labels):
                            if item.isupper(): raw_exp_data[item].append(items[i+1])
                            elif item.islower(): raw_calc_data[item].append(items[i+1])
                            else: _hybs.append(items[i+1])
                    continue
                if line.startswith("EXCH"):
                    items = line.strip().split()[1:]
                    _exchange_list.append(items)

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

    def report(self):
        for _nmrdata in self.nmrdata:
            for i, item in enumerate(_nmrdata.all_exchange_list):
                #temp_list = lambda y: [_nmrdata.label[x] for x in chain.from_iterable(y)]
                print("=="*10+"("+str(i+1)+")"+"=="*10)
                #print(",".join(temp_list(self.all_exchange_list[0])) + "==>" + ",".join(temp_list(item)))
                #print("--"*20)
                temp_exp_data = _nmrdata.exchangeNMR(item)
                self.printNMR(_nmrdata.label, _nmrdata.calc_data, temp_exp_data, _nmrdata.dtype)

if __name__ == "__main__":
    s = StatisticalNMR()
    filename = "../data/aldols.C"
    s.readDataFromFile(filename, '13C')
    filename = "../data/aldols.H"
    s.readDataFromFile(filename, '1H')
    s.report()
