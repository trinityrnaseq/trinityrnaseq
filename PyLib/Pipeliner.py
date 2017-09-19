#!/usr/bin/env python
# encoding: utf-8

#from __future__ import (absolute_import, division,
#                        print_function, unicode_literals)

import os, sys
import logging
import subprocess
import shlex

logger = logging.getLogger(__name__)

def run_cmd(cmd):
    cmd = shlex.split(cmd)
    logger.info("Running: " + " ".join(cmd))
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (output, error) = process.communicate()
    logger.debug(output)
    if process.returncode != 0:
        logger.error(error)
        raise RuntimeError("Error while running command \"" + str(cmd) + "\":\n" + error)

    
class Pipeliner(object):

    _checkpoint_dir = None
    _cmds_list = []

    def __init__(self, checkpoint_dir):

        checkpoint_dir = os.path.abspath(checkpoint_dir)

        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
            
        self._checkpoint_dir = checkpoint_dir
    


    def add_commands(self, cmds_list):

        for cmd in cmds_list:
            self._cmds_list.append(cmd)

    
    def num_cmds(self):
        return len(self._cmds_list)


    def run(self):
        for cmd in self._cmds_list:
            checkpoint_file = os.path.sep.join([self._checkpoint_dir, cmd.get_checkpoint()])
            if os.path.exists(checkpoint_file):
                logger.info("CMD: " + cmd.get_cmd() + " already processed. Skipping.")
            else:
                # execute it.  If it succeeds, make the checkpoint file
                run_cmd(cmd.get_cmd())
                run_cmd("touch {}".format(checkpoint_file))

        # since all commands executed successfully, remove them from the current cmds list
        self._cmds_list = list()
    


class Command(object):

    _cmd = None
    _checkpoint = None

    def __init__(self, cmd, checkpoint):
        self._cmd = cmd
        self._checkpoint = checkpoint

    def get_cmd(self):
        return self._cmd

    def get_checkpoint(self):
        return self._checkpoint

