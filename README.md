# SD4SOLPS
Synthetic diagnostic workflow manager for use with SOLPS models

This repository is the top level layer for managing a workflow for calculating
the outputs of synthetic diagnostics attached to a SOLPS model.

Steps:
1) Load SOLPS into IMAS DD format if not already
2) Load equilibrium (that the SOLPS mesh was based on) into IMAS DD format
3) Make assumptions to extend profiles into the core and far SOL, if needed
4) Run synthetic diagnostic models and record output
