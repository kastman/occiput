
# Import packages and modules 
import numpy as np


# Exceptions 
class BadParameter(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr("Bas parameter: "+str(self.value))


# Print with 3 levels of verbosity 
verbose = 1
def set_verbose_high(): 
    """Print everything - DEBUG mode"""
    global verbose; verbose = 2    
def set_verbose_medium(): 
    """Print runtime information"""
    global verbose; verbose = 1
def set_verbose_low(): 
    """Print only important messages"""
    global verbose; verbose = 0
def set_verbose_no_printing(): 
    """Do not print messages at all"""
    global verbose; verbose = -1        
def get_verbose_level():
    return verbose
def print_high_verbose(msg):
    """Use this for DEBUG Information"""
    if verbose >= 2: 
        print msg
def print_medium_verbose(msg):
    """Use this for messages useful at runtime"""
    if verbose >= 1: 
        print msg
def print_low_verbose(msg):
    """Use this for important messages"""
    if verbose >= 0: 
        print msg


