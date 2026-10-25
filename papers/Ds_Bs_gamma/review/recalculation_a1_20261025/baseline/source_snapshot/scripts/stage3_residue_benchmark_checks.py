"""Compatibility entry point for the corrected Wang validation.

The former overlap-fit exercise compared mixed-current residues with an
obsolete 0.245 GeV value from an early Wang version.  The current arXiv v4
value is 0.345 GeV, and only the pure AA correlator is directly comparable.
"""

from wang_pure_axial_validation import main


if __name__ == "__main__":
    main()
