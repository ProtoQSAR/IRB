# -*- coding: utf-8 -*-

import json
import sys
import os
import pandas as pd
import pickle



def load_json(json_path):
    '''Loads JSON files as a dictionary.
    '''
    if hasattr(sys, '_MEIPASS'):
        json_path = os.path.join(sys._MEIPASS, json_path)
    with open(json_path, 'r', encoding='utf8') as json_file:
            json_dict = json.load(json_file)

    return json_dict




def nou_read_csv(path, **kwargs):
    '''Loads JSON files as a dictionary.
    '''
    if hasattr(sys, '_MEIPASS'):
        path = os.path.join(sys._MEIPASS, path)
            
    
    return pd.read_csv(path, **kwargs)

def nou_open(path, mode, **kwargs):
    if hasattr(sys, '_MEIPASS'):
        path = os.path.join(sys._MEIPASS, path)
    return open(path, mode, **kwargs)

     
def save_json(json_path, diccionario):
    '''Saves JSON files
    '''
    with open(json_path, 'w', encoding='utf8') as json_file:
        json.dump(diccionario, json_file, indent=4)
    
