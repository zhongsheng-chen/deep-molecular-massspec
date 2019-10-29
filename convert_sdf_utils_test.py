#!/usr/bin/python 
# -*- coding: utf-8 -*-  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Author: Zhongsheng Chen 
# Date: 10/22/2019 
# Copyright: Copyright 2019, Beijing University of Chemical Technology 
# License: The MIT License (MIT)
# Email: zhongsheng.chen@outlook.com 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

"""Tests for convert_sdf_util."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import tempfile

import tensorflow as tf
from absl import flags
from absl.testing import absltest
from absl.testing import parameterized

import convert_sdf_utils
import parse_sdf_utils


def _make_test_dir(relative_path):
    return os.path.join(
        flags.FLAGS.test_srcdir, os.path.split(os.path.abspath(__file__))[0], relative_path)


class ConvertSDFUtilsTest(tf.test.TestCase, parameterized.TestCase):

    def setUp(self):
        self.test_data_directory = _make_test_dir('test_dataset/')
        self.test_bad_sdf_name = os.path.join(self.test_data_directory, 'test_mona_vf_npl.sdf')
        self.out_dir = tempfile.mkdtemp(dir=absltest.get_default_test_tmpdir())

    def tearDown(self):
        tf.gfile.DeleteRecursively(self.out_dir)

    def test_convert_to_sdf(self):
        convert_sdf_utils.convert_to_sdf(
            self.test_bad_sdf_name,
            failed_block_file_name='test_mona_vf_npl_failed_blocks.sdf',
            output_dir=self.out_dir)
        path_to_converted_sdf = os.path.join(self.out_dir,
                                             'converted_{}'.format(os.path.split(self.test_bad_sdf_name)[1]))
        mol_list = parse_sdf_utils.get_sdf_to_mol(path_to_converted_sdf)
        expected_mol_list_len = 1
        self.assertEqual(len(mol_list), expected_mol_list_len)


if __name__ == '__main__':
    tf.test.main()
