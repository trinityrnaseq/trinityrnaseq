#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging
import argparse
import collections
import numpy
import time


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

           
class DP_matrix:
    """
    defines the dynamic programming matrix for the node multiple alignments
    """
    
    @staticmethod
    def build_DP_matrix(num_rows, num_cols):
        dp_matrix = list()
        for i in range(0, num_rows+1):
            row = []
            for j in range(0, num_cols+1):
                struct = { 'i' : i,
                           'j' : j,
                           'bt' : -1,
                           'score' : 0,
                           }
                row.append(struct)
            dp_matrix.append(row)
        
        return dp_matrix    


    @staticmethod
    def toString(dp_matrix):
        nrow = len(dp_matrix)
        ncol = len(dp_matrix[0])

        logger.debug("Matrix is {} X {}".format(nrow, ncol))

        ret_text = ""
        for row in dp_matrix:
            str_row = [ str(x['score']) for x in row ]
            ret_text += "\t".join(str_row) + "\n"

        return ret_text

