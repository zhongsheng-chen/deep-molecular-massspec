#!/usr/bin/env python
# -*- coding: utf-8 -*-
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Author: Zhongsheng Chen 
# Date: 10/10/2019 
# Copyright: Copyright 2019, Beijing University of Chemical Technology 
# License: The MIT License (MIT)
# Email: zhongsheng.chen@outlook.com
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""A helper function for processing non-standardized SDF files.

For dataset from Massbank of North America (MoNA), they are SDF-like files but exact SDF files. In the SDF-like file,
some lines line header sections was missing. So, these files can not used be loaded as a standard SDF file using RDKit
Tool. Specifications of a standard SDF file are given at https://en.wikipedia.org/wiki/Chemical_table_file.
Molecule block are loaded and then append 'M  END' and other lines to make sure it can be loaded as sdf.

Example:
        $ python make_standard_sdf.py \
            --path_to_sdf_like_file=/sdf/like/file/path \
            --output_dir_name=/tmp/sdflike_to_sdf \
            --alsologtostderr

"""
import os

import tensorflow as tf
from absl import logging
from absl import flags
from absl import app
from rdkit import Chem

FLAGS = flags.FLAGS
flags.DEFINE_string('path_to_sdf_like_file',
                    'testdata/MONA_2_mend.sdf',
                    'specify a full path of a SDF-like file to be standardized as a SDF file')
flags.DEFINE_string('output_dir_name',
                    '/tmp/sdflike_to_sdf',
                    'specify a directory for the SDF-like file converted to SDF.')


def _make_chunk_from_block(block):
    """ Make a chunk from a molecule block"""
    chunk = []
    for line in block:
        if "V2000" in line:
            chunk.extend(['\r', line])
        elif ">  <NAME>" in line:
            chunk.extend(['M  END', line])
        else:
            chunk.append(line)

    return chunk


def _write_chunk_to_file(save_to_path, chunk, mode='w'):
    """ Write a chunk of molecule block to a file"""
    with tf.gfile.Open(save_to_path, mode) as writer:
        for line in chunk:
            writer.write('%s\n' % line)
        # for line_ in chunk_:
        #     writer.write(line_)


def convert_sdflike_to_sdf(unprocessed_sdf_name, output_dir):
    """ Make a standardized sdf
    To fast grape a block context for each molecule, I will extract lines between four dollar signs($$$$), which is
    a flag of end for each molecule in the SDF format.
    """
    suppl = Chem.SDMolSupplier(unprocessed_sdf_name)

    chunk_list = []
    for index in range(len(suppl)):
        molecule_block = suppl.GetItemText(index).strip().splitlines()

        chunk = _make_chunk_from_block(molecule_block)
        chunk_list.append(chunk)

    _, sdf_name = os.path.split(unprocessed_sdf_name)
    save_sdf_to_path = os.path.join(output_dir, 'converted_' + sdf_name)

    for index, chunk in enumerate(chunk_list):
        if index == 0:
            _write_chunk_to_file(save_sdf_to_path, chunk, mode='w')
        else:
            _write_chunk_to_file(save_sdf_to_path, chunk, mode='a')
    num_molecule = len(suppl)
    logging.warning(('Processing on %s from Massbank of North America (MoNA) finished. '
                     '%d molecules has been converted to a standard SDF saved in the path %s'),
                    sdf_name, num_molecule, save_sdf_to_path)


def main(_):
    tf.gfile.MkDir(FLAGS.output_dir_name)
    convert_sdflike_to_sdf(FLAGS.path_to_sdf_like_file, FLAGS.output_dir_name)


if __name__ == '__main__':
    app.run(main)
