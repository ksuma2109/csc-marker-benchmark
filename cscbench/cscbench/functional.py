"""Functional CSC gates — the label-independent ground truth.

A ``FunctionalGate`` holds a per-gene *functional log2 fold change* (CSC vs.
non-CSC) derived from an assay that separates cells by phenotype rather than
transcript: ALDEFLUOR/ALDH sorting, the CD44+CD24- surface gate, CD133 sorting,
or sphere-vs-adherent self-renewal. It can be built from a raw expression
matrix (with group columns) or from a precomputed fold change (e.g. a DESeq2
table).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence, Optional
import numpy as np
import pandas as pd

_BAD_SYMBOLS = {"nan", "none", "na", ""}


def _clean_index(s: pd.Series) -> pd.Series:
    """Coerce index to clean gene-symbol strings; drop dup/empty labels."""
    s = s.copy()
    s.index = s.index.astype(str)
    s = s[~s.index.duplicated(keep="first")]
    s = s[~s.index.str.lower().isin(_BAD_SYMBOLS)]
    return s.dropna()


@dataclass
class FunctionalGate:
    """A functional CSC assay expressed as per-gene log2FC(CSC / non-CSC).

    Attributes
    ----------
    log2fc : pandas.Series
        Index = gene symbol, value = functional log2 fold change.
    name : str
        Human-readable gate name (e.g. "Prostate ALDH+ (GSE270565)").
    cancer : str, optional
        Cancer type, for grouping (e.g. "Prostate").
    criterion : str, optional
        Functional criterion ("ALDH", "CD44/CD24", "CD133", "sphere", ...).
    """

    log2fc: pd.Series
    name: str
    cancer: Optional[str] = None
    criterion: Optional[str] = None

    def __post_init__(self):
        self.log2fc = _clean_index(pd.Series(self.log2fc))

    # ---- constructors ------------------------------------------------------
    @classmethod
    def from_log2fc(cls, log2fc, name: str, **kw) -> "FunctionalGate":
        """Build from a precomputed per-gene log2 fold change (Series or dict)."""
        return cls(log2fc=pd.Series(log2fc), name=name, **kw)

    @classmethod
    def from_matrix(
        cls,
        expr: pd.DataFrame,
        csc_cols: Sequence[str],
        noncsc_cols: Sequence[str],
        name: str,
        gene_col: Optional[str] = None,
        pseudocount: float = 1.0,
        **kw,
    ) -> "FunctionalGate":
        """Build from an expression matrix and CSC / non-CSC sample columns.

        log2FC = log2(mean(CSC) + c) - log2(mean(non-CSC) + c).

        Parameters
        ----------
        expr : DataFrame
            Rows = genes, columns include the CSC and non-CSC sample columns.
        csc_cols, noncsc_cols : list of str
            Column names for the two groups.
        name : str
            Gate name.
        gene_col : str, optional
            Column holding gene symbols. If None, ``expr.index`` is used.
        pseudocount : float
            Added before the log (default 1.0).
        """
        df = expr
        symbols = df[gene_col] if gene_col is not None else pd.Series(df.index, index=df.index)
        csc = df[list(csc_cols)].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        non = df[list(noncsc_cols)].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        lfc = np.log2(csc + pseudocount) - np.log2(non + pseudocount)
        s = pd.Series(np.asarray(lfc), index=pd.Index(symbols, dtype=str))
        return cls(log2fc=s, name=name, **kw)

    # ---- helpers -----------------------------------------------------------
    @property
    def genes(self) -> set:
        return set(self.log2fc.index)

    def top_quartile_labels(self) -> pd.Series:
        """Binary truth: 1 if a gene is in the top functional quartile."""
        return (self.log2fc > self.log2fc.quantile(0.75)).astype(int)

    def __len__(self) -> int:
        return len(self.log2fc)

    def __repr__(self) -> str:
        tag = f" [{self.cancer}/{self.criterion}]" if self.cancer else ""
        return f"FunctionalGate({self.name!r}{tag}, {len(self)} genes)"
