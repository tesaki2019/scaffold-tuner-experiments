# utils/scaffold_tuner_benchmark.py
#
# Benchmark-specific wrapper for Scaffold-tuner.
# Compatible with Scaffold-tuner v1.0.3 (commit 5bdfc8c).
#
# Uses ONLY the public API of scaffold_tuner.scaffold_intervention:
#   make_murcko_core_with_rlabels, get_r_labels, get_original_rgroups,
#   build_molecule, build_extended_rgroup_smiles, get_fragment_pool,
#   count_features
#
# This module reimplements the generation loop so that candidate_rank is
# assigned BEFORE the descriptor-delta filter, matching the convention used
# by Random, CReM, and LibInvent.  The public propose_structures() is NOT
# modified.
#
# Generation strategies (mirrors v1.0.3 propose_structures behaviour):
#   "tier1" = whole R-group replacement  (remove_* modes)
#   "tier2" = within-R-group extension   (add_* modes via build_extended_rgroup_smiles)

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from scaffold_tuner.scaffold_intervention import (
    make_murcko_core_with_rlabels,
    get_r_labels,
    get_original_rgroups,
    build_molecule,
    build_extended_rgroup_smiles,
    get_fragment_pool,
    count_features,
)

# Descriptor index in count_features() return tuple (v1.0.3):
#   (HBA=0, HBD=1, AR=2, RotB=3, HA=4)
_DESC_INDEX = {"HBA": 0, "HBD": 1, "AR": 2, "ROTB": 3, "HA": 4}

# Canonical mode strings used by v1.0.3
_MODE_FROM_TARGET = {
    ("HBD",  +1): "add_hbd",
    ("HBD",  -1): "remove_hbd",
    ("HBA",  +1): "add_hba",
    ("HBA",  -1): "remove_hba",
    ("AR",   +1): "add_ar",
    ("AR",   -1): "remove_ar",
    ("ROTB", +1): "add_rb",
    ("ROTB", -1): "remove_rb",
}


def _collateral(parent_feats, gen_feats, target_idx):
    """Sum of |delta| for all non-target descriptors (indices 0-3)."""
    return sum(
        abs(gen_feats[i] - parent_feats[i])
        for i in range(4)
        if i != target_idx
    )


def generate_all_candidates(
    parent_smiles: str,
    target=None,
    delta=None,
    mode: str | None = None,
    max_candidates: int = 100,
    random_seed=None,
) -> tuple[list[dict], dict]:
    """Generate Scaffold-tuner candidates for benchmark use (v1.0.3 compatible).

    Returns ALL valid, scaffold-retained, non-duplicate candidates BEFORE
    filtering on descriptor delta.  candidate_rank is assigned at this stage
    so it has the same meaning as in Random, CReM, and LibInvent.

    Parameters
    ----------
    parent_smiles   : SMILES of the parent molecule
    target          : descriptor key — HBD, HBA, AR, RotB (case-insensitive)
    delta           : requested delta (+1 or -1)
    mode            : v1.0.3 mode string (e.g. "add_hba"); overrides target/delta
    max_candidates  : max candidates to return (applied after dedup/scaffold)
    random_seed     : unused, kept for API consistency

    Returns
    -------
    (candidates, stats)
    candidates : list of dicts per accepted candidate, each with:
                   generated_smiles, candidate_rank, generation_strategy,
                   to_fragment, site, actual_delta, collateral_change,
                   scaffold_preserved, parent_features, generated_features
    stats      : {n_examined, n_invalid, n_scaffold_failed, n_duplicate,
                  n_valid_unique}
                 Identity: n_examined = n_invalid + n_scaffold_failed
                                      + n_duplicate + n_valid_unique
    """
    # Resolve mode
    if mode is not None:
        resolved_mode = mode
    elif target is not None and delta is not None:
        key = (target.upper(), int(delta))
        resolved_mode = _MODE_FROM_TARGET.get(key)
        if resolved_mode is None:
            raise ValueError(f"Unsupported (target, delta): {key}")
    else:
        raise ValueError("Provide either mode or (target, delta).")

    direction = +1 if resolved_mode.startswith("add_") else -1
    # Descriptor index from mode name
    desc_map = {
        "hbd": "HBD", "hba": "HBA", "ar": "AR", "rb": "ROTB", "ha": "HA"
    }
    suffix = resolved_mode.split("_", 1)[1]
    target_key = desc_map.get(suffix)
    if target_key is None:
        raise ValueError(f"Cannot parse descriptor from mode: {resolved_mode}")
    target_idx = _DESC_INDEX[target_key]

    # Build scaffold / R-group decomposition
    parent_mol, scaffold_mol, core = make_murcko_core_with_rlabels(parent_smiles)
    parent_canon     = Chem.MolToSmiles(parent_mol)
    core_smiles      = Chem.MolToSmiles(core)
    labels           = get_r_labels(core)
    original_rgroups = get_original_rgroups(parent_smiles, core)
    parent_feats     = count_features(parent_mol)
    parent_scaffold_smi = Chem.MolToSmiles(scaffold_mol)

    # Stats — n_examined incremented BEFORE every build attempt.
    # Identity: n_examined = n_invalid + n_scaffold_failed
    #                       + n_duplicate + n_valid_unique
    n_examined        = 0
    n_invalid         = 0
    n_scaffold_failed = 0
    n_duplicate       = 0

    seen: set[str]        = set()
    all_valid: list[dict] = []

    for label in labels:
        site = f"R{label}"
        pool = get_fragment_pool(resolved_mode, label)

        for frag_smiles in pool:
            n_examined += 1   # count every attempt BEFORE build

            try:
                if direction > 0:
                    # add_* mode: extend the existing substituent
                    # (within-R-group extension = tier2)
                    new_frag = build_extended_rgroup_smiles(
                        original_rgroups[site], frag_smiles, label
                    )
                    if new_frag is None:
                        # Fully substituted site; treat as duplicate/skip
                        n_duplicate += 1
                        continue
                    rgroups   = {**original_rgroups, site: new_frag}
                    strategy  = "tier2"
                    to_frag   = new_frag
                else:
                    # remove_* mode: replace the substituent outright
                    # (whole R-group replacement = tier1)
                    rgroups   = {**original_rgroups, site: frag_smiles}
                    strategy  = "tier1"
                    to_frag   = frag_smiles

                gen_mol   = build_molecule(core_smiles, rgroups)
                gen_canon = Chem.MolToSmiles(gen_mol)
            except Exception:
                n_invalid += 1
                continue

            if gen_canon == parent_canon or gen_canon in seen:
                n_duplicate += 1
                continue

            # Scaffold retention (same criterion as all other methods)
            if not gen_mol.HasSubstructMatch(scaffold_mol):
                n_scaffold_failed += 1
                continue

            # ── candidate accepted ──────────────────────────────────────────
            seen.add(gen_canon)
            gen_feats    = count_features(gen_mol)
            actual_delta = gen_feats[target_idx] - parent_feats[target_idx]
            collateral   = _collateral(parent_feats, gen_feats, target_idx)

            gen_scaffold_mol = MurckoScaffold.GetScaffoldForMol(gen_mol)
            gen_scaffold_smi = (
                Chem.MolToSmiles(gen_scaffold_mol)
                if gen_scaffold_mol is not None else ""
            )

            all_valid.append({
                "generated_smiles":          gen_canon,
                "generation_strategy":       strategy,
                "to_fragment":               to_frag,
                "site":                      site,
                "actual_delta":              actual_delta,
                "collateral_change":         collateral,
                "scaffold_preserved":        True,
                "parent_features":           parent_feats,
                "generated_features":        gen_feats,
                "parent_scaffold_smiles":    parent_scaffold_smi,
                "generated_murcko_scaffold": gen_scaffold_smi,
            })

    # candidate_rank = acceptance order (no collateral/score-based re-sorting)
    ranked = all_valid[:max_candidates]
    for rank, rec in enumerate(ranked, start=1):
        rec["candidate_rank"] = rank

    n_valid_unique = len(all_valid)   # before max_candidates cap

    # Verify identity
    assert n_examined == n_invalid + n_scaffold_failed + n_duplicate + n_valid_unique, (
        f"Stats identity failed: {n_examined} != "
        f"{n_invalid}+{n_scaffold_failed}+{n_duplicate}+{n_valid_unique}"
    )

    return ranked, {
        "n_examined":        n_examined,
        "n_invalid":         n_invalid,
        "n_scaffold_failed": n_scaffold_failed,
        "n_duplicate":       n_duplicate,
        "n_valid_unique":    n_valid_unique,
    }
