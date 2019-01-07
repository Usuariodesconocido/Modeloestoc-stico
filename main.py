import configparser
import logging
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, './src')
from BaseModel import BaseModel

def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                logging.info('skip: %s' % option)
        except:
            logging.info('exception on %s!' % option)
            dict1[option] = None
    return dict1


def get_general_conf(name):
    Config = configparser.ConfigParser()
    Config.read('./conf/config.txt')
    myprior = {}
    for sec in Config.sections():
        if sec == name:
            myprior = ConfigSectionMap(sec, Config)
    return myprior

if __name__ == '__main__':
    generalconf = get_general_conf('Mconf')
    autorregresiveconf = get_general_conf('ARconf')

    logging.basicConfig(
        format='%(asctime)s %(levelname)s:%(message)s',
        level=logging.DEBUG)
    logging.info('Modelo estocastico')

    base = BaseModel(Mconf=generalconf, ARconf=autorregresiveconf)
    M = base.calculate_matrix_M(Mconf=generalconf)
    VV= base.value_function(Mconf=generalconf)
    F=base.value_function(Mconf=generalconf)
    plt.show(F)