#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

logger = logging.getLogger(__name__)


def get_Trinity_gene_name(transcript_name):
    """
    extracts the gene name from the Trinity identifier as the prefix
    """

    (gene_name, count) = re.subn("_i\d+$", "", transcript_name)
    if count != 1:
        errmsg = "Error, couldn't extract gene_id from transcript_id: {}".format(transcript_name)
        logger.critical(errmsg)
        raise RuntimeError(errmsg)

    return gene_name


