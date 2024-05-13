#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------------

import pandas as pd
import numpy as np

def updateMeta(config):
    _sanity_checks(config)
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', dtype=str)
    #print(metadata)    
    #LEN: For the colummns below return a dictionary where key is col name
    #sans prefixes, and value = ([list of samples w/ 1], [list of samples w/2])
    for col in ['mergeRNA', 'mergeATAC', 'multiome']:
        #Get all columns of that type
        tmp = filter(lambda x: x.startswith(col+"_"), metadata.columns)
        for c in tmp:
            #LEN: UPDATE- only one group allowed
            samples_lst1 = list(metadata.loc[metadata[c] == '1'].index.values)
            #samples_lst2 = list(metadata.loc[metadata[c] == '2'].index.values)
            col_name = c[len(col)+1:]
            if col in config:
                #config[col][col_name] = (samples_lst1, samples_lst2)
                config[col][col_name] = samples_lst1
            else:
                #config[col] = {col_name: (samples_lst1, samples_lst2)}
                config[col] = {col_name: samples_lst1}
    batch_1 = list(metadata.loc[metadata['batch'] == '1'].index.values)
    batch_2 = list(metadata.loc[metadata['batch'] == '2'].index.values)
    config['batch'] = (batch_1, batch_2)

    assay_RNA = list(metadata.loc[metadata['assay'] == 'RNA'].index.values)
    assay_ATAC = list(metadata.loc[metadata['assay'] == 'ATAC'].index.values)
    config['assay'] = (assay_RNA, assay_ATAC)
    config["ordered_sample_list"] = metadata.index
    return config


def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', ' ':'_', '/':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()

