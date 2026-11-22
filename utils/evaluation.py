# utils/evaluation.py
# Success / strict-success judgement and Success@k calculation.
# Import these functions rather than re-implementing per notebook.

import pandas as pd
from configs.benchmark_config import TASKS, EVAL_DESCRIPTORS


def calculate_descriptor_delta(parent_desc, candidate_desc):
    """Return per-descriptor delta dict: candidate - parent."""
    return {d: candidate_desc[d] - parent_desc[d] for d in EVAL_DESCRIPTORS}


def is_success(parent_desc, candidate_desc, task_name):
    """True if target descriptor changed by exactly the requested delta."""
    task = TASKS[task_name]
    actual = candidate_desc[task["descriptor"]] - parent_desc[task["descriptor"]]
    return actual == task["delta"]


def is_strict_success(parent_desc, candidate_desc, task_name):
    """True if target descriptor changed by requested delta AND all others unchanged."""
    task = TASKS[task_name]
    target_desc = task["descriptor"]
    if candidate_desc[target_desc] - parent_desc[target_desc] != task["delta"]:
        return False
    return all(
        candidate_desc[d] == parent_desc[d]
        for d in EVAL_DESCRIPTORS if d != target_desc
    )


def success_at_k(df, k, all_parent_ids, metric="success"):
    """Compute Success@k per (parent_id, task), with 0-candidate pairs filled as False.

    Parameters
    ----------
    df            : DataFrame with columns [parent_id, task, candidate_rank, success/strict_success]
    k             : int   candidate budget
    all_parent_ids: iterable of all parent_ids in the benchmark
    metric        : "success" or "strict_success"

    Returns
    -------
    Series indexed by (parent_id, task) with bool values.
    Missing (parent_id, task) pairs → False (denominator is never inflated).
    """
    from configs.benchmark_config import TASK_NAMES

    # Complete index: all parent × task combinations
    full_index = pd.MultiIndex.from_product(
        [sorted(all_parent_ids), TASK_NAMES],
        names=["parent_id", "task"],
    )

    result = (
        df[df["candidate_rank"] <= k]
        .groupby(["parent_id", "task"])[metric]
        .any()
        .reindex(full_index, fill_value=False)
    )
    return result


def success_rate_at_k(df, k, all_parent_ids, metric="success"):
    """Overall Success@k: mean over all (parent_id, task) pairs."""
    return success_at_k(df, k, all_parent_ids, metric=metric).mean()


def success_rate_by_task_at_k(df, k, all_parent_ids, metric="success"):
    """Task-wise Success@k: mean over parents, per task."""
    s = success_at_k(df, k, all_parent_ids, metric=metric)
    return s.groupby("task").mean()
