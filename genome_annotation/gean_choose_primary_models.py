#!/usr/bin/env python

import argparse
import pysam
import sys
from collections import defaultdict

"""gean_choose_primary_models.py:
Given an input gene annotation file, choose the primary model (i.e. mRNA with
the longest ORF) and drop all alternate models ."""

