#!/usr/bin/env python


# Script to create transposed pedfile from VCF file

# The MIT License (MIT)
#
# Copyright (c) 2015 Beaulieu-Saucier Pharmacogenomics Centre
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from __future__ import print_function

import os
import re
import sys
import gzip
import shlex
import logging
import argparse
import unittest
import collections
from shutil import copyfile
from collections import defaultdict

import pandas as pd


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Louis-Philippe Lemieux Perreault. "
                 "All rights reserved.")
__license__ = "MIT"
__credits__ = ["Louis-Philippe Lemieux Perreault", ]
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


# The version of the script
__version__ = "0.4"


# The testing mode
_TESTING_MODE = False


# Configuring the log
logger = logging.getLogger("vcf2tped")


def main(args=None):
    """The main function."""
    # The parser object
    desc = """Convert a VCF to a TPED (version {}).""".format(__version__)
    parser = argparse.ArgumentParser(description=desc)

    # Running the script
    try:
        # Getting and checking the options
        args = parse_args(parser, args)
        check_args(args)

        if not _TESTING_MODE:
            sh = logging.StreamHandler()
            sh.setFormatter(logging.Formatter(
                fmt="[%(asctime)s %(levelname)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            ))
            logger.addHandler(sh)
            logger.setLevel(logging.INFO)

            logger.info("This is {script_name} version {version}".format(
                script_name=os.path.basename(sys.argv[0]),
                version=__version__,
            ))
            logger.info("Arguments: {}".format(
                " ".join(shlex.quote(part) for part in sys.argv[1:]),
            ))

        # Reading the ped file
        sample_info = read_ped(args.ped)

        # Converting
        convert_vcf(args.vcf, sample_info, args.out, args.skip_haploid)

    except KeyboardInterrupt:
        logger.info("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    except ProgramError as e:
        logger.error(e.message)
        if not _TESTING_MODE:
            parser.error(e.message)
        sys.exit(1)


def read_ped(i_filename):
    """Reads the PED file to gather sample information.

    :param i_filename: the name of the input file
    :type i_filename: str

    :returns: sample information

    """
    # Checking the number of columns of the file
    with open(i_filename, "r") as i_file:
        if len(i_file.readline().split("\t")) != 6:
            msg = "{}: wrong number of column".format(i_filename)
            raise ProgramError(msg)

    data = pd.read_csv(i_filename, sep="\t",
                       names=["FID", "IID", "Father", "Mother", "Gender",
                              "Status"])

    return data.set_index("IID", drop=False, verify_integrity=True)


def convert_vcf(i_filename, sample_info, o_prefix, skip_haploid):
    """Reads a VCF and creates a TPED/TFAM from it.

    :param i_filename: the name of the VCF file (might be gzip).
    :param sample_info: information about the samples.
    :param o_prefix: the prefix of the output files.
    :param skip_haploid: whether to check haploid genotypes or not.

    :type i_filename: string
    :type sample_info: pandas.DataFrame
    :type o_prefix: string
    :type skip_haploid: bool

    """
    # Some regular expression
    single_point_re = re.compile("^\.$")
    allele_split_re = re.compile("[/|]")

    # The open function to use
    open_f = open
    if re.search("\.gz$", i_filename):
        open_f = gzip.open

    # Reading the file
    with open_f(i_filename, 'rb') as i_file:
        line = i_file.readline().decode()

        # We read until the header
        while re.search("^##", line) is not None:
            line = i_file.readline().decode()

        # This should be the header
        if re.search("^##CHROM", line) is not None:
            msg = "{}: no header".format(i_filename)
            raise ProgramError(msg)

        # Creating the header
        row = line.rstrip("\r\n").split("\t")
        header = {name:  i for i, name in enumerate(row)}

        # Checking some names
        for name in ["#CHROM", "POS", "ID", "REF", "ALT", "FORMAT"]:
            if name not in header:
                msg = "{}: no column named {}".format(i_filename, name)
                raise ProgramError(msg)

        # Printing the TFAMs
        tfam_names = ["{}.snv.2_alleles.tfam".format(o_prefix),
                      "{}.snv.n_alleles.tfam".format(o_prefix),
                      "{}.indel.2_alleles.tfam".format(o_prefix),
                      "{}.indel.n_alleles.tfam".format(o_prefix)]
        sample_info = print_same_tfams(tfam_names, sample_info,
                                       row[header["FORMAT"]+1:])

        # Those positions have been seen
        seen_pos = defaultdict(int)

        # The output files
        tped_snv_2 = None
        tped_snv_n = None
        tped_indel_2 = None
        tped_indel_n = None
        snv_ref = None
        try:
            tped_snv_2 = open("{}.snv.2_alleles.tped".format(o_prefix), 'w')
            tped_snv_n = open("{}.snv.n_alleles.tped".format(o_prefix), 'w')
            tped_indel_2 = open(
                "{}.indel.2_alleles.tped".format(o_prefix),
                "w",
            )
            tped_indel_n = open(
                "{}.indel.n_alleles.tped".format(o_prefix),
                "w",
            )
            snv_ref = open("{}.snv.ref".format(o_prefix), "w")
            indel_ref = open("{}.indel.ref".format(o_prefix), "w")
        except IOError:
            msg = "couldn't write output files"
            raise ProgramError(msg)

        # Reading the rest of the data
        for line in i_file:
            row = line.decode().rstrip("\r\n").split("\t")

            # Getting the information
            chrom = encode_chr(row[header["#CHROM"]])
            pos = row[header["POS"]]
            name = row[header["ID"]]
            alleles = [row[header["REF"]]] + row[header["ALT"]].split(",")
            g_format = row[header["FORMAT"]].split(":")
            g_format = {name: i for i, name in enumerate(g_format)}
            genotypes = [
                "." if single_point_re.match(i)
                else i.split(":")[g_format["GT"]]
                for i in row[header["FORMAT"]+1:]
            ]

            # Getting rid of the "." (so that it becomes "./.")
            genotypes = [single_point_re.sub("./.", i) for i in genotypes]

            # Is this an INDEL?
            indel = False
            for allele in alleles:
                if len(allele) > 1:
                    indel = True
                    break

            # The output file to choose from (default is SNV and bi-allelic)
            o_file = tped_snv_2
            o_ref_file = snv_ref
            if len(alleles) > 2:
                o_file = tped_snv_n

            # Constructing the genotype map
            g_map = {str(i): a for i, a in enumerate(alleles)}
            if indel:
                # This is a new genotype map, only for INDELs
                g_map = {str(i): str(i+1) for i in range(len(alleles))}

                # These are the required files
                o_ref_file = indel_ref
                o_file = tped_indel_2
                if len(alleles) > 2:
                    o_file = tped_indel_n

            # Adding the unknown genotype
            g_map["."] = "0"

            # Checking if the position have a name
            if name == ".":
                name = "{}:{}".format(chrom, pos)

            # We increase the number of time we saw this name, and check if we
            # saw it more than once
            seen_pos[name] += 1
            if seen_pos[name] > 1:
                name = "{}-{}".format(name, seen_pos[name])

            # The first part of the TPED line
            first_part = [chrom, name, "0", pos]

            genotypes = [allele_split_re.split(i) for i in genotypes]
            genotypes = [
                recode_genotype(g, g_map, chrom, pos, sample_info.iloc[i, 0],
                                sample_info.iloc[i, 4], skip_haploid)
                for i, g in enumerate(genotypes)
            ]
            print("\t".join(first_part + genotypes), file=o_file)

            # Saving the alleles
            print("\t".join([chrom, pos, name, alleles[0],
                             ",".join(alleles[1:])]),
                  file=o_ref_file)

        # Closing the output files
        tped_snv_2.close()
        tped_snv_n.close()
        tped_indel_2.close()
        tped_indel_n.close()
        snv_ref.close()
        indel_ref.close()


def recode_genotype(genotype, encoding, chromosome, position, sample, gender,
                    skip_haploid):
    """Encodes the genotypes.

    :param genotype: the genotypes (list of alleles)
    :param encoding: the allele encoding
    :param chromosome: the chromosome on which the marker is
    :param position: the position of the marker
    :param sample: the ID of the sample
    :param gender: the gender of the sample
    :param skip_haploid: whether to check or not the haploid genotypes

    :type genotype: list of str
    :type encoding: map
    :type chromosome: str
    :type position: str
    :type sample: str
    :type gender: int
    :type skip_haploid: bool

    :returns: a string with the two alleles separated by a space.

    """
    if len(genotype) == 1:
        # This is an haploid genotype
        if gender == 1 and (chromosome == "23" or chromosome == "24"):
            return " ".join(encoding[genotype[0]] * 2)

        elif skip_haploid:
            logger.warning(
                "chr{chrom}:{pos}: {sample} (gender {gender}): haploid call "
                "set as homozygous".format(
                    chrom=chromosome, pos=position, sample=sample,
                    gender=gender,
                ),
            )
            return " ".join(encoding[genotype[0]] * 2)

        else:
            logger.warning(
                "chr{chrom}:{pos}: {sample} (gender {gender}): haploid call "
                "set as no call".format(
                    chrom=chromosome, pos=position, sample=sample,
                    gender=gender,
                ),
            )
            return " ".join([encoding["."]] * 2)

    return "{} {}".format(encoding[genotype[0]], encoding[genotype[1]])


def print_same_tfams(file_names, sample_info, sample_names):
    """Print the same TFAM in different files.

    :param file_names: list of output file names
    :param sample_info: information about the samples
    :param sample_names: the list of samples

    :type file_names: list of string
    :type sample_info: pandas.DataFrame
    :type sample_names: list of string

    :returns: the final TFAM

    """
    # Checking the samples
    if sorted(sample_names) != sorted(list(sample_info.IID)):
        msg = "samples in PED are not the same in the VCF"
        raise ProgramError(msg)

    # Reordering the sample information
    sample_info = sample_info.loc[sample_names, :]

    try:
        # Printing only the first TFAM
        sample_info.to_csv(file_names[0], sep="\t", header=False,
                           index=False)

        # Copying the first TFAM into the others
        for name in file_names[1:]:
            copyfile(file_names[0], name)
    except IOError:
        msg = "couldn't write output TFAMs"
        raise ProgramError(msg)

    return sample_info


def encode_chr(chrom):
    """Encodes a chromosome.

    :param chrom: the chromosome to encode
    :type chrom: string

    Remove the leading "chr" if necessary, and encodes sex chromosomes and
    other to their numeric counterpart.

        - X    X chromosome                    -> 23
        - Y    Y chromosome                    -> 24
        - XY   Pseudo-autosomal region of X    -> 25
        - MT   Mitochondrial                   -> 26

    .. doctest::

        >>> [encode_chr(str(i)) for i in range(0, 11)]
        ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
        >>> [encode_chr("chr{}".format((i))) for i in range(11, 21)]
        ["11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
        >>> [encode_chr(str(i)) for i in range(21, 27)]
        ["21", "22", "23", "24", "25", "26"]
        >>> [encode_chr(i) for i in ["X", "Y", "XY", "MT"]]
        ["23", "24", "25", "26"]
        >>> [encode_chr("chr{}".format(i) for i in range(21, 27)]
        ["23", "24", "25", "26"]
        >>> encode_chr("27")
        Traceback (most recent call last):
            ...
        ProgramError: 27: not a valid chromosome
        >>> encode_chr("0")
        Traceback (most recent call last):
            ...
        ProgramError: 0: not a valid chromosome
        >>> encode_chr("chrXX")
        Traceback (most recent call last):
            ...
        ProgramError: chrXX: not a valid chromosome

    """
    # Removing the leading "chr"
    chrom = re.sub("^chr", "", chrom, flags=re.I)

    # Checking for special chromosome
    if re.match("^x$", chrom, flags=re.I):
        return "23"
    if re.match("^y$", chrom, flags=re.I):
        return "24"
    if re.match("^xy$", chrom, flags=re.I):
        return "25"
    if re.match("^mt$", chrom, flags=re.I):
        return "26"

    try:
        numerical_chrom = int(chrom)
        if numerical_chrom >= 1 and numerical_chrom <= 26:
            return chrom
    except ValueError:
        pass

    msg = "{}: not a valid chromosome".format(chrom)
    raise ProgramError(msg)


def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    if not os.path.isfile(args.vcf):
        msg = "{}: no such file".format(args.vcf)
        raise ProgramError(msg)

    if not ((re.search("\.vcf$", args.vcf)
            or re.search("\.vcf\.gz$", args.vcf))):
        msg = "{}: not a vcf file".format(args.vcf)
        raise ProgramError(msg)

    if not os.path.isfile(args.ped):
        msg = "{}: no such file".format(args.ped)
        raise ProgramError(msg)

    return True


def parse_args(parser, args=None):
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    parser.add_argument(
        "-v", "--version", action="version",
        version="{} version {}".format(os.path.basename(sys.argv[0]),
                                       __version__),
    )

    # The input file
    group = parser.add_argument_group("Input Files")
    group.add_argument(
        "--vcf", type=str, metavar="FILE", required=True,
        help="The VCF file.",
    )
    group.add_argument(
        "--ped", type=str, metavar="FILE", required=True,
        help="The PED file.",
    )

    # Program options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--skip-haploid-check", action="store_true", dest="skip_haploid",
        help="If used, no check will be performed for haploid genotypes. They "
             "will be converted to homozygous of the same allele.",
    )

    # The output file
    group = parser.add_argument_group("Output Files")
    group.add_argument(
        "-o", "--out", type=str, metavar="STR", default="vcf_to_tped",
        help="The suffix of the output file. [Default: %(default)s]",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


class TestVCF2TPED(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        global _TESTING_MODE
        _TESTING_MODE = True

    def setUp(self):
        """Setting up the test cases."""
        # Setting up the parser information
        self.prefix = os.path.join("test", "vcf2tped_test")
        args = [
            "--vcf", os.path.join("test", "input.vcf"),
            "--ped", os.path.join("test", "input.ped"),
            "--out", self.prefix,
        ]

        self.suffixes = [".indel.2_alleles.tfam", ".indel.2_alleles.tped",
                         ".indel.n_alleles.tfam", ".indel.n_alleles.tped",
                         ".indel.ref", ".snv.2_alleles.tfam",
                         ".snv.2_alleles.tped", ".snv.n_alleles.tfam",
                         ".snv.n_alleles.tped", ".snv.ref"]

        # Executing the script
        main(args=args)

        # Now the test for the gender
        self.prefix_gender = os.path.join("test", "vcf2tped_gender_test")
        args = [
            "--vcf", os.path.join("test", "input.gender.vcf"),
            "--ped", os.path.join("test", "input.gender.ped"),
            "--out", self.prefix_gender,
        ]

        # Executing the script
        with self._assertLogs() as logs:
            main(args=args)
        self.logs = logs.output

    def tearDown(self):
        """Deletes the output files."""
        import glob

        # Deleting files for the main test
        for filename in glob.glob(self.prefix + ".*"):
            os.remove(filename)

        # Deleting files for the gender test
        for filename in glob.glob(self.prefix_gender + ".*"):
            os.remove(filename)

    def test_output_file_present(self):
        """Check if all output files are present."""
        for filename in [self.prefix + suffix for suffix in self.suffixes]:
            self.assertTrue(os.path.isfile(filename))

    def test_input_file_missing(self):
        """Tests the script raises an error if an input file is missing."""
        args = [
            "--vcf", os.path.join("test", "wrong_input.vcf"),
            "--ped", os.path.join("test", "input.ped"),
            "--out", self.prefix,
        ]
        with self.assertRaises(SystemExit) as cm:
            with self._assertLogs(level=logging.ERROR) as logs:
                main(args=args)

        # Checking the exception
        self.assertNotEqual(cm.exception.code, 0)

        # Checking the error log
        self.assertEqual(
            ["ERROR:vcf2tped:test/wrong_input.vcf: no such file"],
            logs.output,
        )

        args = [
            "--vcf", os.path.join("test", "input.vcf"),
            "--ped", os.path.join("test", "wrong_input.ped"),
            "--out", self.prefix,
        ]
        with self.assertRaises(SystemExit) as cm:
            with self._assertLogs(level=logging.ERROR) as logs:
                main(args=args)

        # Checking the exception
        self.assertNotEqual(cm.exception.code, 0)

        # Checking the error log
        self.assertEqual(
            ["ERROR:vcf2tped:test/wrong_input.ped: no such file"],
            logs.output,
        )

    def test_missing_samples_in_pedfile(self):
        """Tests the script raises an error sample missing in pedfile."""
        args = [
            "--vcf", os.path.join("test", "input.vcf"),
            "--ped", os.path.join("test", "missing_samples.ped"),
            "--out", self.prefix,
        ]
        with self.assertRaises(SystemExit) as cm:
            with self._assertLogs(level=logging.ERROR) as logs:
                main(args=args)

        # Checking the exception
        self.assertNotEqual(cm.exception.code, 0)

        # Checking the error log
        self.assertEqual(
            ["ERROR:vcf2tped:samples in PED are not the same in the VCF"],
            logs.output,
        )

    def test_not_a_vcf_file(self):
        """Tests when the input file is not a vcf file."""
        args = [
            "--vcf", os.path.join("test", "input.notvcf"),
            "--ped", os.path.join("test", "input.ped"),
            "--out", self.prefix,
        ]

        # Creating the vcf file
        with open(os.path.join("test", "input.notvcf"), "w"):
            pass

        with self.assertRaises(SystemExit) as cm:
            with self._assertLogs(level=logging.ERROR) as logs:
                main(args=args)

        # Deleting the VCF file
        os.remove(os.path.join("test", "input.notvcf"))

        # Checking the exception
        self.assertNotEqual(cm.exception.code, 0)

        # Checking the error log
        self.assertEqual(
            ["ERROR:vcf2tped:test/input.notvcf: not a vcf file"],
            logs.output,
        )

    def test_invalid_chromosomes(self):
        """Tests an invalid chromosome."""
        # Invalid int chromosome
        with self.assertRaises(ProgramError) as cm:
            encode_chr("0")
        self.assertEqual("0: not a valid chromosome", cm.exception.message)

        # Invalid str chromosome
        with self.assertRaises(ProgramError) as cm:
            encode_chr("XYZ")
        self.assertEqual("XYZ: not a valid chromosome", cm.exception.message)

    def test_xy_chromosome(self):
        """Tests the XY chromosome."""
        self.assertEqual("25", encode_chr("XY"))

    def test_tfams(self):
        output = (
            "HG00096\tHG00096\t0\t0\t0\t-9\n"
            "HG00097\tHG00097\t0\t0\t0\t-9\n"
            "HG00099\tHG00099\t0\t0\t0\t-9\n"
            "HG00100\tHG00100\t0\t0\t0\t-9\n"
        )
        for filename in [self.prefix + suffix for suffix in self.suffixes]:
            if not filename.endswith(".tfam"):
                continue

            content = None
            with open(filename, "r") as i_file:
                content = i_file.read()

            self.assertEqual(output, content)

    def test_snv_2_alleles(self):
        """Tests the output for SNV with 2 alleles."""
        output = (
            "24\trs139377059\t0\t16050678\tC C\tC T\tC T\tT T\n"
            "22\trs6518357\t0\t16051107\tC C\tC A\tC A\tC C\n"
            "22\t22:16051453\t0\t16051453\tG G\t0 0\tC G\tG C\n"
            "22\t22:16051453-2\t0\t16051453\tA A\tC A\tC A\tA A\n"
            "23\trs79725552\t0\t16051480\tT T\tT C\tT C\tT T\n"
            "23\trs79725552-2\t0\t16051480\tT T\tT C\tT C\tT T\n"
            "26\trs79725552-3\t0\t16051480\tT T\tT C\tT C\tT T\n"
        )

        content = None
        with open(self.prefix + ".snv.2_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_snv_n_alleles(self):
        """Tests the output for SNV with more than 2 alleles."""
        output = (
            "22\trs188945759\t0\t16050984\tC A\tC C\tC G\tG A\n"
        )

        content = None
        with open(self.prefix + ".snv.n_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_snv_ref(self):
        """Checks the content of the SNV ref file."""
        output = (
            "24\t16050678\trs139377059\tC\tT\n"
            "22\t16050984\trs188945759\tC\tG,A\n"
            "22\t16051107\trs6518357\tC\tA\n"
            "22\t16051453\t22:16051453\tG\tC\n"
            "22\t16051453\t22:16051453-2\tA\tC\n"
            "23\t16051480\trs79725552\tT\tC\n"
            "23\t16051480\trs79725552-2\tT\tC\n"
            "26\t16051480\trs79725552-3\tT\tC\n"
        )

        content = None
        with open(self.prefix + ".snv.ref", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_indel_2_alleles(self):
        """Tests the output for INDEL with 2 alleles."""
        output = (
            "22\trs149201999\t0\t16050408\t1 1\t1 2\t1 2\t1 1\n"
            "22\trs146752890\t0\t16050612\t1 2\t1 2\t0 0\t1 1\n"
            "22\t22:16051453-3\t0\t16051453\t1 1\t1 1\t1 1\t0 0\n"
            "1\trs1234\t0\t16051453\t1 1\t1 1\t1 1\t0 0\n"
        )

        content = None
        with open(self.prefix + ".indel.2_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_indel_n_alleles(self):
        """Tests the output for INDEL with n alleles."""
        output = (
            "22\trs62224609\t0\t16051249\t1 1\t2 1\t2 1\t2 3\n"
        )

        content = None
        with open(self.prefix + ".indel.n_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_indel_ref(self):
        """Checks the content of the INDEL ref file."""
        output = (
            "22\t16050408\trs149201999\tT\tCC\n"
            "22\t16050612\trs146752890\tCC\tC\n"
            "22\t16051249\trs62224609\tT\tTC,TCC\n"
            "22\t16051453\t22:16051453-3\tC\tCA\n"
            "1\t16051453\trs1234\tC\tCA\n"
        )

        content = None
        with open(self.prefix + ".indel.ref", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_gender(self):
        """Tests effect of gender on genotypes."""
        output = (
            "22\trs149201\t0\t16050408\t0 0\t0 0\t0 0\n"
            "23\trs149202\t0\t16050408\t0 0\t0 0\tT T\n"
            "24\trs149203\t0\t16050408\t0 0\t0 0\tT T\n"
            "24\trs149204\t0\t16050408\tT T\tT C\tT C\n"
            "23\trs149205\t0\t16050408\tT T\tT C\tT C\n"
        )

        content = None
        with open(self.prefix_gender + ".snv.2_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_gender_skip_haploid_check(self):
        """Tests effect of gender on genotypes when skipping haploid check."""
        args = [
            "--vcf", os.path.join("test", "input.gender.vcf"),
            "--ped", os.path.join("test", "input.gender.ped"),
            "--out", self.prefix_gender,
            "--skip-haploid-check"
        ]
        with self._assertLogs(level=logging.WARNING) as logs:
            main(args=args)

        # Checking the log
        self.assertEqual(
            ["WARNING:vcf2tped:chr22:16050408: HG00096 (gender 2): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr22:16050408: HG00097 (gender 0): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr22:16050408: HG00099 (gender 1): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr23:16050408: HG00096 (gender 2): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr23:16050408: HG00097 (gender 0): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr24:16050408: HG00096 (gender 2): haploid "
             "call set as homozygous",
             "WARNING:vcf2tped:chr24:16050408: HG00097 (gender 0): haploid "
             "call set as homozygous"],
            logs.output,
        )

        # Checking the output
        output = (
            "22\trs149201\t0\t16050408\tT T\tC C\tT T\n"
            "23\trs149202\t0\t16050408\tT T\tC C\tT T\n"
            "24\trs149203\t0\t16050408\tT T\tT T\tT T\n"
            "24\trs149204\t0\t16050408\tT T\tT C\tT C\n"
            "23\trs149205\t0\t16050408\tT T\tT C\tT C\n"
        )

        content = None
        with open(self.prefix_gender + ".snv.2_alleles.tped", "r") as i_file:
            content = i_file.read()

        self.assertEqual(output, content)

    def test_logs(self):
        """Checks the logs."""
        self.assertEqual(
            ["WARNING:vcf2tped:chr22:16050408: HG00096 (gender 2): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr22:16050408: HG00097 (gender 0): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr22:16050408: HG00099 (gender 1): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr23:16050408: HG00096 (gender 2): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr23:16050408: HG00097 (gender 0): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr24:16050408: HG00096 (gender 2): haploid "
             "call set as no call",
             "WARNING:vcf2tped:chr24:16050408: HG00097 (gender 0): haploid "
             "call set as no call"],
            self.logs,
        )

    def _assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)


class BaseTestCaseContext_Compatibility:

    def __init__(self, test_case):
        self.test_case = test_case

    def _raiseFailure(self, standardMsg):
        msg = self.test_case._formatMessage(self.msg, standardMsg)
        raise self.test_case.failureException(msg)


_LoggingWatcher = collections.namedtuple("_LoggingWatcher",
                                         ["records", "output"])


class CapturingHandler_Compatibility(logging.Handler):
    """
    A logging handler capturing all (raw and formatted) logging output.
    """

    def __init__(self):
        logging.Handler.__init__(self)
        self.watcher = _LoggingWatcher([], [])

    def flush(self):
        pass

    def emit(self, record):
        self.watcher.records.append(record)
        msg = self.format(record)
        self.watcher.output.append(msg)


class AssertLogsContext_Compatibility(BaseTestCaseContext_Compatibility):
    """A context manager used to implement TestCase.assertLogs()."""

    LOGGING_FORMAT = "%(levelname)s:%(name)s:%(message)s"

    def __init__(self, test_case, logger_name, level):
        BaseTestCaseContext_Compatibility.__init__(self, test_case)
        self.logger_name = logger_name
        if level:
            # Python < 3.4 logging doesn't have a _nameToLevel dictionary
            nameToLevel = {
                'CRITICAL': logging.CRITICAL,
                'ERROR': logging.ERROR,
                'WARN': logging.WARNING,
                'WARNING': logging.WARNING,
                'INFO': logging.INFO,
                'DEBUG': logging.DEBUG,
                'NOTSET': logging.NOTSET,
            }
            self.level = nameToLevel.get(level, level)
        else:
            self.level = logging.INFO
        self.msg = None

    def __enter__(self):
        if isinstance(self.logger_name, logging.Logger):
            logger = self.logger = self.logger_name
        else:
            logger = self.logger = logging.getLogger(self.logger_name)
        formatter = logging.Formatter(self.LOGGING_FORMAT)
        handler = CapturingHandler_Compatibility()
        handler.setFormatter(formatter)
        self.watcher = handler.watcher
        self.old_handlers = logger.handlers[:]
        self.old_level = logger.level
        self.old_propagate = logger.propagate
        logger.handlers = [handler]
        logger.setLevel(self.level)
        logger.propagate = False
        return handler.watcher

    def __exit__(self, exc_type, exc_value, tb):
        self.logger.handlers = self.old_handlers
        self.logger.propagate = self.old_propagate
        self.logger.setLevel(self.old_level)
        if exc_type is not None:
            # let unexpected exceptions pass through
            return False
        if len(self.watcher.records) == 0:
            self._raiseFailure(
                "no logs of level {} or higher triggered on {}"
                .format(logging.getLevelName(self.level), self.logger.name))


# Calling the main, if necessary
if __name__ == "__main__":
    main()
