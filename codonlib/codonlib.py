"""Main module."""
from collections import defaultdict
from functools import cache
from itertools import product
from operator import itemgetter
from typing import List

import numpy as np

from Bio.Data.CodonTable import unambiguous_dna_by_id


class CodonDegeneracy:
    def __init__(self, table_id: int = 1):
        self.table_id = table_id
        self.codontable_atlas = unambiguous_dna_by_id[self.table_id]
        self.codon2aa = {}
        self.aa2codons = defaultdict(list)
        self.codon_table: np.char.array = None
        self.nt_table: np.char.array = None
        self.aa_table: np.char.array = None
        self.__codon_aa_mappings()

    def __codon_aa_mappings(self):
        for codon, aa in self.codontable_atlas.forward_table.items():
            self.codon2aa[codon] = aa
            self.aa2codons[aa].append(codon)
        for codon in self.codontable_atlas.stop_codons:
            self.codon2aa[codon] = "*"
            self.aa2codons["*"].append(codon)

    def __create_tables(self):
        codon_list = []
        aa_list = []
        for i, col_nt in enumerate(["T", "C", "A", "G"]):
            for j, wobble_nt in enumerate(["T", "C", "A", "G"]):
                for k, row_nt in enumerate(["T", "C", "A", "G"]):
                    codon = col_nt + row_nt + wobble_nt
                    aa = self.codon2aa[codon]
                    codon_list.append(codon)
                    aa_list.append(aa)
        self.codon_table = np.char.array(codon_list, dtype=str).reshape((16, 4))
        self.nt_table = self.codon_table.view("U1").reshape((16, 4, -1))
        self.aa_table = np.char.array(aa_list).reshape((16, 4))

    @cache
    def __(self, nt):
        print(nt)
        codon_list = [x + y + z for x in nt[0] for y in nt[1] for z in nt[2]]
        return set(itemgetter(*codon_list)(self.codon2aa))

    @cache
    def __get_aa_possibilities(self, codons: list) -> set:
        """
        Get the possible combinations from a given list of possibilities.

        Parameters
        ----------
        nt_combo : list
            List of codons.

        Examples
        --------
        >>>__locknkey(['CCG', 'AAG'])
        {'T', 'K', 'P', 'Q'}

        Returns
        -------
        set
            amino acid symbol set
        """
        nt1, nt2, nt3 = frozenset(), frozenset(), frozenset()
        for codon in codons:
            nt1 |= frozenset(codon[0])
            nt2 |= frozenset(codon[1])
            nt3 |= frozenset(codon[2])
        return self.__combinations((nt1, nt2, nt3))

    def off_targets(self, aa_list: List[str]) -> set:
        """
        Get the off-target amino acids for a given list of amino acids.

        Parameters
        ----------
        aa_list : list
            List of amino acid codons.

        Returns
        -------
        set
            Set of off-target amino acids.
        """
        on_target_aa = set(aa_list)
        off_target_aa: set = set()
        off_target_best = [0] * 42
        all_codons_per_aa = [self.aa2codons[aa] for aa in aa_list]
        for one_codon_selected_per_aa in product(*all_codons_per_aa):
            aa_possibilities = self.__get_aa_possibilities(one_codon_selected_per_aa)
            off_target_aa = aa_possibilities - on_target_aa
            if len(off_target_aa) < len(off_target_best):
                off_target_best = off_target_aa
        # needs to return best combos and list of all equals of each combos so i know the "real" set of codons to return
        # in real life we want all the combons not just 1 for each aa
        return off_target_best
