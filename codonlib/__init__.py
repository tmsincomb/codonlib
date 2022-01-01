"""Top-level package for codonlib."""
from .codonlib import CodonDegeneracy


codon_degeneracy = CodonDegeneracy()
off_targets = codon_degeneracy.off_targets

__author__ = """Troy Sincomb"""
__email__ = "tsincomb@gmail.com"
__version__ = "0.1.0"
__all__ = ["off_targets", "codon_degeneracy"]
