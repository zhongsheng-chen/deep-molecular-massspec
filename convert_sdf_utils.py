#!/usr/bin/env python
# -*- coding: utf-8 -*-
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Author: Zhongsheng Chen
# Date: 10/10/2019
# Copyright: Copyright 2019, Beijing University of Chemical Technology
# License: The MIT License (MIT)
# Email: zhongsheng.chen@outlook.com
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

r"""A helper function for processing SDF-like dataset files from MoNA.

For dataset from Massbank of North America (MoNA), they are SDF-like files but exact SDF files. In the SDF-like file,
some lines in header sections was missing. So, these files can not used be loaded as a standard SDF file using RDKit
Tool. Specifications of a standard SDF file are given at https://en.wikipedia.org/wiki/Chemical_table_file.
Molecule block are loaded and then append 'M  END' and other lines to make sure it can be loaded as sdf.

Example:
        $ python convert_sdf_utils.py \
            --path_to_bad_sdf=/sdf/like/file/path \
            --failed_block_file_name=/save/failed/block/to/file \
            --output_dir=/save/path/to/converted/sdf \
            --alsologtostderr
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

import numpy as np
import tensorflow as tf
from absl import app
from absl import flags
from absl import logging
from openbabel import pybel
from rdkit import Chem

FLAGS = flags.FLAGS
flags.DEFINE_string('path_to_bad_sdf',
                    'test_dataset/test_mona_vf_npl.sdf',
                    'specify a relative path of a SDF-like file to convert as SDF file')
flags.DEFINE_string('failed_block_file_name',
                    '',
                    'specify a file to store failed molecule blocks')
flags.DEFINE_string('output_dir',
                    '',
                    'specify a directory for SDF files converted.')

# required tags in sdf files.
SDF_TAG_MASS_SPEC_PEAKS = 'MASS SPECTRAL PEAKS'
SDF_TAG_INCHIKEY = 'INCHIKEY'
SDF_TAG_INCHI = 'INCHI'
SDF_TAG_NAME = 'NAME'
SDF_TAG_MOLECULE_MASS = 'EXACT MASS'

expected_props = [SDF_TAG_MASS_SPEC_PEAKS,
                  SDF_TAG_INCHIKEY,
                  SDF_TAG_INCHI,
                  SDF_TAG_NAME,
                  SDF_TAG_MOLECULE_MASS]


def _make_mol_block_from_string(mol_str_in_lines):
    """ Make a molecule block from a string read by Chem.SDMolSupplier

    A valid molecule header should have three separated lines: title line, program line and  counts line.
    The missing lines in the original dataset files are filled up with new blank lines (''). To restored molecules
    from a SDF data file, it needs to include 'M  END' as molecule reading end flag. Due to google's
    Neural Electron−Ionization Mass Spectrometry (NEIMS) demanding InChIKey, I extracted InChI from
    the comment section in SDF data files, and then transfer them to the corresponding InChiKey,
    and write them back to molecule blocks in SDF data files.

    Args:
        mol_str_in_lines: List of string for molecule descriptions.

    Returns:
        A list of string for molecule block.
    """

    block = []
    for line in mol_str_in_lines:
        if "V2000" in line:  # insert a required lines in molecule header block.
            block.extend(['', line])
        elif ">  <NAME>" in line:  # mark M  END as an token for a valid molecule.
            block.extend(['M  END', line])
        elif "$$$$" in line:  # add up two properties, including INCHIKEY and INCHI before ending a molecule block.
            inchikey = ''
            inchi = _get_prop_value_from_mol_block(block, SDF_TAG_INCHI)
            if inchi:
                inchikey = Chem.inchi.InchiToInchiKey('InChI=' + inchi)
                # logging.warning('InChI %s has a InChiKey %s', inchi, inchikey)
                # uncomment for logging InChi and its InChiKey
            block.extend(['>  <%s>' % 'INCHI', inchi, ''])
            block.extend(['>  <%s>' % SDF_TAG_INCHIKEY, inchikey, '', line])
        else:
            block.append(line)
    return block


def _make_mol_block_rational(mol_block):
    """ Make molecule blocks passed parse successfully.

    Note that although some missing lines and a molecule reading flag (M END) has been filled up,
    for some unknown reasons, considerable molecule blocks can not be successfully recognized as expected.
    I found that OpenBabel can calibrate those. The reason this phenomenon may
    lie in wrong values in their atom blocks values, making molecule blocks corrupted.

    Args:
        mol_block: A raw mol_block (list of string )

    Returns
        A molecule block whose atom block is calibrated, if molecule block is successfully loaded as sdf,
        otherwise, None.
    """

    mol_str = '\n'.join(mol_block)
    try:
        mol_obj = pybel.readstring('sdf', mol_str)
    except IOError:
        mol_obj = None

    return mol_obj.write('sdf').splitlines() if mol_obj is not None else None


def _has_prop_on_mol_block(block, prop_key):
    """ Check if a specific property exits in the molecule block.

    All properties sored in SDF data files include but not limited to
    SDF_TAG_MASS_SPEC_PEAKS, SDF_TAG_INCHIKEY,
    SDF_TAG_INCHI, SDF_TAG_NAME,
    SDF_TAG_MOLECULE_MASS

    Args:
        block: A molecule block (string).
        prop_key: A key of a property of molecule.

    Returns:
        True if has the specific property prop_key, otherwise, False.

    Raises:
        ValueError if prop_key does not match any element in the expected properties.
    """

    if prop_key not in expected_props:
        raise ValueError('%s is not a supported property type.', prop_key)
    has_prop = False
    for line in block:
        if line.strip() == ('>  <%s>' % prop_key):
            has_prop = True
    return has_prop


def _get_prop_value_from_mol_block(block, prop_key):
    """ Get the corresponding value for a specific property of any molecule.

    Assert that each InChI for a molecule is stored in the first line of COMMENT section of a SDF data file which
    start with 'InChI=InChI'.

    Args:
        block: A list of string to specify a molecule block.
        prop_key: A string to specify key of a property of molecules.

    Returns:
        value of the specific property prop_key. if the specific property of the molecule block does not find,
        return ''.
    """

    if _has_prop_on_mol_block(block, prop_key) or prop_key == SDF_TAG_INCHI:
        prop_key_full_name = ('>  <%s>' % prop_key.upper())
        comment_key_full_name = '>  <COMMENT>'
        prop_value = ''
        if prop_key_full_name == '>  <INCHI>':
            prop_key_full_name = '>  <COMMENT>'
        for ind, line in enumerate(block):

            if prop_key_full_name == comment_key_full_name and prop_key_full_name == line.strip():
                if block[ind + 1].startswith('InChI='):
                    prop_value = block[ind + 1].strip('InChI=')
                    break
            elif prop_key_full_name == line.strip():
                prop_value = block[ind + 1].strip('')
                break
        return prop_value


def _write_mol_block_to_file(save_to_path, mol_block_list):
    """ Write a molecule block to a file"""

    def __mol_block_writer(mode):
        with tf.gfile.Open(save_to_path, mode) as writer:
            for line in mol_block:
                writer.write('%s\n' % line)

    for index, mol_block in enumerate(mol_block_list):
        if index == 0:
            __mol_block_writer(mode='w')
        else:
            __mol_block_writer(mode='a')


def _check_mol_block_has_all_prop(mol_block):
    """ Check if all the considered properties exist in the molecule block."""
    prop_check_status = []
    for mol_prop_key in expected_props:
        prop_check_status.append(_has_prop_on_mol_block(mol_block, mol_prop_key))

    # I will check if each block of every molecule has all the tags:
    #  'MASS SPECTRAL PEAKS', 'INCHIKEY', 'INCHI', NAME', 'EXACT MASS'
    return np.all(prop_check_status)


def convert_to_sdf(path_to_bad_sdf, failed_block_file_name=None, output_dir=None):
    """ Make a sdf-like data file converted to sdf.

    Note that to form a proper molecule block, some missing lines in header of original molecule string read from
    sdf-like data files are filled up. Moreover, InChI will be taken to convert to its InChIKey. Both
    InChI and InChIKey are taken as properties and then append them leading to four dollar signs ($$$$)
    in each molecule string.

    Args:
        path_to_bad_sdf: Path to load a sdf-like data file.
        failed_block_file_name: Specify a file to failed molecule blocks, if set. By default (None), failed molecule
            blocks are skipped, no failed molecule block written to file.
        output_dir: Directory to the converted sdf file if set.  By default (None), the directory of the sdf-like data
            file will be set as output directory.
    """

    suppl = Chem.SDMolSupplier(path_to_bad_sdf)

    mol_block_list = []
    logging.warning('Converting started ...')
    for index in range(len(suppl)):
        # all string leading to each four dollar signs($$$$) will take as molecule-related string.
        mol_str = suppl.GetItemText(index).strip().splitlines()
        mol_block = _make_mol_block_rational(_make_mol_block_from_string(mol_str))

        if mol_block is not None:
            mol_block_list.append(mol_block)

    num_failed = 0
    valid_mol_block_list = []
    failed_mol_block_list = []
    for mol_block in mol_block_list:
        if _check_mol_block_has_all_prop(mol_block):
            valid_mol_block_list.append(mol_block)
        else:
            num_failed += 1
            failed_mol_block_list.append(mol_block)

    out_dir, sdf_name = os.path.split(path_to_bad_sdf)
    if output_dir is '':
        output_dir = os.path.abspath(out_dir)
    output_dir = os.path.abspath(output_dir)
    save_valid_mol_block_to_path = os.path.join(output_dir, ('converted_%s' % sdf_name))
    save_failed_mol_block_to_path = os.path.join(output_dir, failed_block_file_name)

    if valid_mol_block_list:
        _write_mol_block_to_file(save_valid_mol_block_to_path, valid_mol_block_list)
    if failed_block_file_name and failed_mol_block_list:
        _write_mol_block_to_file(save_failed_mol_block_to_path, failed_mol_block_list)
    num_valid_mol_block = len(valid_mol_block_list)
    num_failed_mol_block = len(failed_mol_block_list)
    logging.warning(('Processing on %s from Massbank of North America (MoNA) finished. '
                     'Except for %d failed molecule blocks, totally, '
                     '%d molecules have been converted to a read-friendly SDF saved in the path %s'),
                    sdf_name, num_failed_mol_block, num_valid_mol_block, save_valid_mol_block_to_path)


def main(_):
    tf.gfile.MkDir(FLAGS.output_dir)
    convert_to_sdf(FLAGS.path_to_bad_sdf, FLAGS.failed_block_file_name, FLAGS.output_dir)


if __name__ == '__main__':
    app.run(main)
