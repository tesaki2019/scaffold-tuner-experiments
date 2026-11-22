# utils/scaffold_tuner_benchmark.py
#
# Benchmark-specific wrapper for Scaffold-tuner v7.
#
# Public API (`propose_structures`) returns only success candidates
# (actual_delta == requested_delta), so candidate_rank there counts only
# successes.  For a fair benchmark we need candidate_rank to reflect the
# position in the full pool of valid, scaffold-retained, non-duplicate
# candidates — BEFORE the descriptor delta filter — exactly as Random,
# CReM, and LibInvent work.
#
# This module reimplements the generation loop from propose_structures,
# keeping all candidates that pass:
#   RDKit build  →  scaffold retention  →  deduplication
# and assigns candidate_rank at that stage.  Descriptor success/strict-success
# is then evaluated by make_candidate_record() as usual.
#
# The public propose_structures() is NOT modified.

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from scaffold_tuner.scaffold_intervention import (
    make_murcko_core_with_rlabels,
    get_r_labels,
    get_original_rgroups,
    _build_molecule,
    _adapt_to_site,
    _rgroup_sub_substitutions,
    normalize_fragment_smiles,
    count_features,
    _collateral_change,
    R_GROUP_SMILES,
    _DESC_INDEX,
    _MODE_MAP,
    _SUPPORTED,
)


def generate_all_candidates(
    parent_smiles: str,
    target: str,
    delta: int,
    mode: str | None = None,
    max_candidates: int = 100,
    random_seed: int | None = None,
) -> tuple[list[dict], dict]:
    """Generate Scaffold-tuner candidates for benchmark use.

    Returns ALL valid, scaffold-retained, non-duplicate candidates
    WITHOUT filtering on descriptor delta.  candidate_rank is assigned
    at this stage so it matches the convention used by Random, CReM,
    and LibInvent.

    Each record in the returned list has:
        generated_smiles, candidate_rank, generation_strategy (tier1/tier2),
        to_fragment, site, actual_delta, collateral_change,
        scaffold_preserved, parent_features, generated_features

    Parameters
    ----------
    parent_smiles   : SMILES of the parent molecule
    target          : descriptor key — HBD, HBA, AR, RotB
    delta           : requested delta (+1 or -1)
    mode            : legacy mode string (overrides target/delta)
    max_candidates  : max candidates to return (applied after dedup/scaffold)
    random_seed     : unused here but kept for API consistency

    Returns
    -------
    (candidates, stats)
    candidates : list of dicts (all valid+scaffold-retained, ranked)
    stats      : {n_examined, n_invalid, n_scaffold_failed, n_duplicate,
                  n_valid_unique}
                 where n_examined counts every pool substitution attempted
    """
    # Resolve mode → target/delta
    if mode is not None and target is None:
        target, delta = _MODE_MAP[mode]

    # Pre-compute pool canonical set for Tier 1/2 classification
    pool_canon = {
        Chem.MolToSmiles(Chem.MolFromSmiles(s))
        for s in R_GROUP_SMILES
        if Chem.MolFromSmiles(s) is not None
    }

    target_upper = target.upper()
    target_idx = _DESC_INDEX.get(target_upper)
    if target_idx is None:
        raise ValueError(f"Unknown target descriptor: {target}")

    # Build scaffold / R-group decomposition
    parent_mol, scaffold_mol, core = make_murcko_core_with_rlabels(parent_smiles)
    parent_canon     = Chem.MolToSmiles(parent_mol)
    core_smiles      = Chem.MolToSmiles(core)
    labels           = get_r_labels(core)
    original_rgroups = get_original_rgroups(parent_smiles, core)
    parent_feats     = count_features(parent_mol)
    parent_scaffold_smi = Chem.MolToSmiles(scaffold_mol)

    # Stats counters.
    # n_examined = every generation attempt (incremented before try/except),
    # so that:  n_examined = n_invalid + n_scaffold_failed + n_duplicate + n_valid_unique
    n_examined        = 0
    n_invalid         = 0   # _build_molecule exception
    n_scaffold_failed = 0
    n_duplicate       = 0

    seen: set[str]        = set()
    all_valid: list[dict] = []

    for label in labels:
        site           = f"R{label}"
        current_rgroup = original_rgroups.get(site, "")

        tier1 = [
            (_adapt_to_site(smi, label), normalize_fragment_smiles(smi), "tier1")
            for smi in R_GROUP_SMILES
        ]
        tier2_rgroups = (
            _rgroup_sub_substitutions(current_rgroup, R_GROUP_SMILES)
            if current_rgroup else []
        )
        tier2 = [
            (s, normalize_fragment_smiles(s), "tier2")
            for s in tier2_rgroups
        ]

        for adapted, to_frag_label, strategy in tier1 + tier2:
            n_examined += 1   # count every attempt before build

            try:
                rgroups   = {**original_rgroups, site: adapted}
                gen_mol   = _build_molecule(core_smiles, rgroups)
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
            collateral   = _collateral_change(parent_feats, gen_feats, target_idx)

            gen_scaffold_mol = MurckoScaffold.GetScaffoldForMol(gen_mol)
            gen_scaffold_smi = (
                Chem.MolToSmiles(gen_scaffold_mol)
                if gen_scaffold_mol is not None else ""
            )

            all_valid.append({
                "generated_smiles":          gen_canon,
                "generation_strategy":       strategy,
                "to_fragment":               to_frag_label,
                "site":                      site,
                "actual_delta":              actual_delta,
                "collateral_change":         collateral,   # saved but NOT used for ranking
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

    # Verify identity: n_examined = n_invalid + n_scaffold_failed + n_duplicate + n_valid_unique
    assert n_examined == n_invalid + n_scaffold_failed + n_duplicate + n_valid_unique, (
        f"Stats mismatch: {n_examined} != "
        f"{n_invalid}+{n_scaffold_failed}+{n_duplicate}+{n_valid_unique}"
    )

    stats = {
        "n_examined":        n_examined,
        "n_invalid":         n_invalid,
        "n_scaffold_failed": n_scaffold_failed,
        "n_duplicate":       n_duplicate,
        "n_valid_unique":    n_valid_unique,
    }

    return ranked, stats
