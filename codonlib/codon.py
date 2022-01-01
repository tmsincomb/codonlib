# from Bio.Data.CodonTable import ambiguous_dna_by_id


class Codon:
    def __init__(self, codon) -> None:
        self.codon = codon

    def __str__(self) -> None:
        return self.codon
