import crdb

def print_column_names(tab):
    for icol, col_name in enumerate(tab.dtype.fields):
        print('%2i' % icol, col_name)

def dump_datafile(quantity, energyType, expName, subExpName, filename, combo_level=0):
    filename = 'lake/' + filename
    print(f'search for {quantity} as a function of {energyType} measured by {expName}')
    
    tab = crdb.query(quantity, energy_type=energyType, combo_level=combo_level, energy_convert_level=0, exp_dates=expName)
 
    subExpNames = set(tab["sub_exp"])
    print('number of datasets found : ', len(subExpNames))
    print(subExpNames)

    adsCodes = set(tab["ads"])
    print(adsCodes)

    items = [i for i in range(len(tab["sub_exp"])) if tab["sub_exp"][i] == subExpName]
    print('number of data : ', len(items))
    assert(len(items) > 0)

    print(f'dump on {filename}')
    with open(filename, 'w') as f:
        f.write(f'#source: CRDB\n')
        f.write(f'#Quantity: {quantity}\n')
        f.write(f'#EnergyType: {energyType}\n')
        f.write(f'#Experiment: {expName}\n')
        f.write(f'#ads: {tab["ads"][items[0]]}\n')
        f.write(f'#E_lo - E_up - y - errSta_lo - errSta_up - errSys_lo - errSys_up\n')
        for eBin, value, errSta, errSys in zip(tab["e_bin"][items], tab["value"][items], tab["err_sta"][items], tab["err_sys"][items]):
            f.write(f'{eBin[0]:10.5e} {eBin[1]:10.5e} {value:10.5e} {errSta[0]:10.5e} {errSta[1]:10.5e} {errSys[0]:10.5e} {errSys[1]:10.5e}\n')
    f.close()
    print('')

if __name__== '__main__':
    dump_datafile('e-', 'R', 'AMS02', 'AMS02 (2011/05-2018/05)', 'AMS-02_e-_rigidity.txt')
    dump_datafile('e+', 'R', 'AMS02', 'AMS02 (2011/05-2018/05)', 'AMS-02_e+_rigidity.txt')
    dump_datafile('H', 'R', 'AMS02', 'AMS02 (2011/05-2018/05)', 'AMS-02_H_rigidity.txt')
    dump_datafile('1H-bar', 'R', 'AMS02', 'AMS02 (2011/05-2018/05)', 'AMS-02_ap_rigidity.txt')
    dump_datafile('1H/e-', 'R', 'AMS02', 'AMS02 (2011/05-2015/05)', 'AMS-02_H_over_e-_rigidity.txt')
    dump_datafile('e-+e+', 'R', 'AMS02', 'AMS02 (2011/05-2018/05)', 'AMS-02_e-e+_rigidity.txt')
    dump_datafile('e-+e+', 'ETOT', 'DAMPE', 'DAMPE (2015/12-2017/06)', 'DAMPE_e-e+_energy.txt')
    dump_datafile('e-+e+', 'ETOT', 'CALET', 'CALET (2015/10-2017/06)', 'CALET_e-e+_energy.txt')
    dump_datafile('e-+e+', 'ETOT', 'FERMI', 'Fermi-LAT-HE (2008/08-2015/06)', 'FERMI_e-e+_energy.txt')

    dump_datafile('e+/e-+e+', 'EK', 'AMS02', 'AMS02 (2011/05-2017/11)', 'AMS-02_pf_energy.txt')
    dump_datafile('e+/e-+e+', 'EK', 'PAMELA', 'PAMELA (2006/07-2009/12)', 'PAMELA_pf_energy.txt')
    dump_datafile('e+/e-+e+', 'EK', 'FERMI', 'Fermi-LAT  (2008/06-2011/04)', 'FERMI_pf_energy.txt')
