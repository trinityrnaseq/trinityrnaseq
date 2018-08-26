#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys
import logging
import subprocess
import shlex
import shutil
import time

logger = logging.getLogger(__name__)


def run_cmd(cmd, ignore_error=False):

    logger.info("Running: " + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        if ignore_error:
            return
        else:
            raise e


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
            # check it's a proper Command object
            if not isinstance(cmd, Command):
                errmsg = "Pipeliner::add_commmands - Error, cmd {} is not a Command object".format(cmd)
                logger.critical(errmsg)
                raise(errmsg)
            
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
                run_cmd(cmd.get_cmd(), cmd.get_ignore_error_setting())
                run_cmd("touch {}".format(checkpoint_file))

        # since all commands executed successfully, remove them from the current cmds list
        self._cmds_list = list()
    


class Command(object):

    def __init__(self, cmd, checkpoint, ignore_error=False):
        self._cmd = cmd
        self._checkpoint = checkpoint
        self._ignore_error = ignore_error

    def get_cmd(self):
        return self._cmd

    def get_checkpoint(self):
        return self._checkpoint

    def get_ignore_error_setting(self):
        return self._ignore_error
 



if __name__ == '__main__':

    checkpoint_dir = "/tmp/checkpoints_dir." + str(time.time())
    
    pipeliner = Pipeliner(checkpoint_dir)

    pipeliner.add_commands([Command("echo hello!", "hello.ok")])

    pipeliner.add_commands([Command("echo done testing pipeliner", "test.ok")])
    
    pipeliner.run()

    shutil.rmtree(checkpoint_dir)

    sys.exit(0)
    
