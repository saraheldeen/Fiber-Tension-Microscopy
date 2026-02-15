#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class FiberTensionParams:
    """
    Parameters for the odd-mode fluctuation model.

    Attributes
    ----------
    lp_um:
        Persistence length ℓp in micrometers [µm].
    kBT_pN_um:
        Thermal energy kBT in [pN·µm] (≈ 0.004114 at ~298 K).
    n_max:
        Largest odd mode index used in the series. Sum uses n = 1,3,...,n_max.
        Must be an odd integer >= 1.
    """
    lp_um: float = 73.89952755
    kBT_pN_um: float = 0.004114
    n_max: int = 999

    def __post_init__(self) -> None:
        if self.lp_um <= 0:
            raise ValueError("lp_um must be positive.")
        if self.kBT_pN_um <= 0:
            raise ValueError("kBT_pN_um must be positive.")
        if self.n_max < 1:
            raise ValueError("n_max must be >= 1.")
        # force odd
        if self.n_max % 2 == 0:
            object.__setattr__(self, "n_max", self.n_max - 1)


class FiberTensionModel:
    """
    Compute filament tension from measured transverse fluctuations (MSD)
    using an odd-mode sum model.

    Model
    -----
    msd = (2/π^4) * (L^3/ℓp) * S(φ)

    with
    S(φ) = Σ_{n odd} 1 / (n^4 + φ n^2)

    and tension (pN):
    τ = (φ π² κ) / L²,  κ = ℓp kBT

    Notes
    -----
    - φ is a dimensionless tension parameter.
    - The series is truncated at n_max (odd modes).
    - φ is solved by bracketed bisection because S(φ) is monotone decreasing.
    """

    def __init__(self, params: FiberTensionParams) -> None:
        self.params = params
        self.pi = float(np.pi)

        # odd modes: 1,3,5,...,n_max
        self.odd_ns = np.arange(1, self.params.n_max + 1, 2, dtype=np.float64)

        # S(0) = sum_{odd} 1/n^4 = π^4/96
        self.S0_limit = (self.pi**4) / 96.0

    def S_of_phi(self, phi: float) -> float:
        """
        S(φ) = Σ_{odd n} 1 / (n^4 + φ n^2).
        """
        n2 = self.odd_ns**2
        return float(np.sum(1.0 / (n2**2 + phi * n2)))

    def phi_from_msd(self, L_um: float, msd_um2: float) -> float:
        """
        Invert MSD -> φ by solving:

            msd = (2/π^4) * (L^3/ℓp) * S(φ)

        Equivalent target:
            S(φ) = (π^4 ℓp msd) / (2 L^3)

        Returns
        -------
        phi:
            Dimensionless φ. Returns np.nan for invalid/unphysical input,
            np.inf for extremely small MSD (effectively φ -> ∞).
        """
        if not (L_um > 0) or not (msd_um2 >= 0):
            return float("nan")

        lp = self.params.lp_um
        S_target = (self.pi**4) * lp * msd_um2 / (2.0 * (L_um**3))

        # If target exceeds S(0), no nonnegative tension can satisfy the model.
        if S_target > self.S0_limit * 1.0000001:
            return float("nan")

        def f(phi: float) -> float:
            return self.S_of_phi(phi) - S_target

        lo, hi = 0.0, 1.0
        flo = f(lo)
        if not np.isfinite(flo):
            return float("nan")

        fhi = f(hi)
        tries = 0
        while fhi > 0.0 and tries < 60:
            hi *= 2.0
            fhi = f(hi)
            tries += 1

        if fhi > 0.0:
            # Could not bracket root: msd -> 0 => φ -> ∞
            return float("inf")

        # Bisection
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = f(mid)

            if abs(fm) < 1e-13:
                return mid

            if fm > 0.0:
                lo = mid
            else:
                hi = mid

            if (hi - lo) < 1e-12 * (1.0 + lo + hi):
                return 0.5 * (lo + hi)

        return 0.5 * (lo + hi)

    def tau_from_phi(self, phi: float, L_um: float) -> float:
        """
        Convert φ to tension τ [pN]:

            τ = (φ π² κ)/L²,  κ = ℓp kBT

        Returns np.nan if inputs are invalid.
        """
        if not np.isfinite(phi) or not (L_um > 0):
            return float("nan")

        kappa = self.params.lp_um * self.params.kBT_pN_um  # [pN·µm²]
        return float((phi * (self.pi**2) * kappa) / (L_um**2))

    def msd_max_tau0(self, L_um: float) -> float:
        """
        Zero-tension MSD ceiling used for diagnostics:

            MSD_max(τ=0) = L^3/(48 ℓp)
        """
        return float((L_um**3) / (48.0 * self.params.lp_um))

    def compute_one(
        self, L_um: float, msd_um2: float
    ) -> Tuple[float, bool, float, float, str]:
        """
        Compute diagnostic values and tension for one (L, MSD) pair.

        Returns
        -------
        msd_max, feasible, phi, tau, note
        """
        msd_max = self.msd_max_tau0(L_um)
        feasible = bool(msd_um2 <= msd_max * 1.0000001)

        phi = self.phi_from_msd(L_um, msd_um2)
        tau = self.tau_from_phi(phi, L_um)

        note = ""
        if not feasible:
            note = "MSD too large for any nonnegative τ (check units/data)."
        elif not np.isfinite(phi):
            note = "Numerical bracket failed (tiny MSD ⇒ huge φ)."

        # Convert inf to nan for cleaner tabular output
        if not np.isfinite(phi):
            phi = float("nan")
        if not np.isfinite(tau):
            tau = float("nan")

        return msd_max, feasible, phi, tau, note

    def compute_dataframe(
        self,
        df: pd.DataFrame,
        length_col: str = "Fiber Length (um)",
        msd_col: str = "Fluctuations (um^2)",
    ) -> pd.DataFrame:
        """
        Compute tension for each row of a dataframe.

        Parameters
        ----------
        df:
            Input dataframe containing length and MSD columns.
        length_col:
            Column name for fiber length in µm.
        msd_col:
            Column name for fluctuations (MSD) in µm².
        """
        if length_col not in df.columns or msd_col not in df.columns:
            raise KeyError(f"Expected columns '{length_col}' and '{msd_col}'.")

        out = df.copy()

        results = []
        for L, msd in zip(out[length_col].astype(float), out[msd_col].astype(float)):
            results.append(self.compute_one(float(L), float(msd)))

        res = pd.DataFrame(
            results,
            columns=[
                "MSD_max at tau=0 (um^2)",
                "Feasible",
                "phi (dimensionless)",
                "Tension (pN)",
                "Note",
            ],
        )
        return pd.concat([out, res], axis=1)


def load_excel_columns(
    xlsx_path: str,
    length_col: Optional[str] = None,
    msd_col: Optional[str] = None,
    fallback_indices: Tuple[int, int] = (1, 3),
    engine: str = "openpyxl",
) -> pd.DataFrame:
    """
    Load an Excel file and extract the relevant columns.

    Parameters
    ----------
    xlsx_path:
        Path to the Excel file.
    length_col, msd_col:
        If provided, use these column names.
        If not provided, use fallback_indices (0-based) like the original script.
    fallback_indices:
        0-based column indices to pick if names aren't provided. Default (1,3).
    """
    df = pd.read_excel(xlsx_path, engine=engine)

    if length_col and msd_col and length_col in df.columns and msd_col in df.columns:
        out = df[[length_col, msd_col]].copy()
        out.columns = ["Fiber Length (um)", "Fluctuations (um^2)"]
        return out

    # Match your original behavior: take columns 2 and 4 (0-based 1 and 3)
    out = df.iloc[:, [fallback_indices[0], fallback_indices[1]]].copy()
    out.columns = ["Fiber Length (um)", "Fluctuations (um^2)"]
    return out

