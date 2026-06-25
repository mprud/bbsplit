"""
bbsplit.parallel
================
Parallel helpers for CPU-bound chemistry work (RDKit parsing, substructure
matching, fragmentation, recombination). Uses processes rather than threads
because RDKit work is CPU-bound and holds the GIL.

How parallelism works here
--------------------------
* A ProcessPoolExecutor spawns independent OS processes. Each process runs on
  whatever core the OS schedules it to, so work is spread across all available
  cores on a multi-core machine. (On a single-core host you still get separate
  processes, just time-sliced onto the one core.)
* Worker functions are module-level (picklable) so the pool works under both
  the 'fork' (Linux) and 'spawn' (Windows/macOS) start methods.
* The RuleSet is shipped to workers as a lightweight, picklable list of dicts,
  because compiled RDKit SMARTS patterns do not pickle reliably. Each worker
  compiles its own SplitEngine once and caches it for the rest of the batch.

`workers` semantics (all public functions)
    workers = 1   -> run serially in the current process (no pool overhead)
    workers = N>1 -> use N worker processes
    workers = 0   -> use os.cpu_count() (all cores)
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from typing import Sequence

from .core import (
    RuleSet, Rule, SplitEngine, Disconnection,
    combine_two,
)


# --------------------------------------------------------------------------- #
#  Picklable rule specs
# --------------------------------------------------------------------------- #
def ruleset_to_spec(ruleset: RuleSet) -> list[dict]:
    """Convert a RuleSet into a plain, picklable list of dicts."""
    return [
        {"name": r.name, "smarts": r.smarts, "priority": r.priority,
         "description": r.description, "enabled": r.enabled}
        for r in ruleset.rules
    ]


def spec_to_ruleset(spec: list[dict]) -> RuleSet:
    return RuleSet([Rule(**d) for d in spec])


def resolve_workers(workers: int) -> int:
    """0 -> all cores; otherwise the requested number (>=1).
    Note: we intentionally do NOT clamp to cpu_count(), so a user may request
    more workers than cores if they wish."""
    if workers is None or workers == 0:
        return os.cpu_count() or 1
    return max(1, int(workers))


# --------------------------------------------------------------------------- #
#  Per-process engine cache (initialised once per worker)
# --------------------------------------------------------------------------- #
_WORKER_ENGINE: SplitEngine | None = None
_WORKER_SIG: tuple | None = None


def _worker_signature(rule_spec: list[dict], max_bonds: int) -> tuple:
    # A content-based signature so the engine is rebuilt only when rules change.
    return (max_bonds, tuple((d["name"], d["smarts"], d["priority"], d["enabled"])
                             for d in rule_spec))


def _ensure_engine(rule_spec: list[dict], max_bonds: int) -> SplitEngine:
    global _WORKER_ENGINE, _WORKER_SIG
    sig = _worker_signature(rule_spec, max_bonds)
    if _WORKER_ENGINE is None or _WORKER_SIG != sig:
        _WORKER_ENGINE = SplitEngine(spec_to_ruleset(rule_spec), max_bonds=max_bonds)
        _WORKER_SIG = sig
    return _WORKER_ENGINE


# --------------------------------------------------------------------------- #
#  Splitting (disconnection) in parallel
# --------------------------------------------------------------------------- #
def _split_one(args) -> tuple[str, list]:
    mol_id, smiles, rule_spec, max_bonds = args
    engine = _ensure_engine(rule_spec, max_bonds)
    return mol_id, engine.disconnect(smiles)


def split_molecules(molecules: Sequence[tuple[str, str]],
                    ruleset: RuleSet,
                    max_bonds: int = 2,
                    workers: int = 1,
                    chunksize: int = 8) -> dict[str, list[Disconnection]]:
    """
    Disconnect many molecules, optionally in parallel across processes.

    molecules : sequence of (mol_id, smiles)
    Returns   : {mol_id: [Disconnection, ...]}
    """
    n = resolve_workers(workers)
    if n == 1 or len(molecules) <= 1:
        engine = SplitEngine(ruleset, max_bonds=max_bonds)
        return {mid: engine.disconnect(smi) for mid, smi in molecules}

    rule_spec = ruleset_to_spec(ruleset)
    tasks = [(mid, smi, rule_spec, max_bonds) for mid, smi in molecules]
    results: dict[str, list[Disconnection]] = {}
    with ProcessPoolExecutor(max_workers=n) as ex:
        for mid, discs in ex.map(_split_one, tasks, chunksize=chunksize):
            results[mid] = discs
    return results


# --------------------------------------------------------------------------- #
#  Enumeration (recombination) in parallel
# --------------------------------------------------------------------------- #
def _combine_pair(pair) -> dict | None:
    i, a, j, b = pair
    prod = combine_two(a, b)
    if prod is None:
        return None
    return {"a_idx": i, "b_idx": j, "blockA": a, "blockB": b, "product": prod}


def enumerate_blocks_parallel(pool_a: Sequence[str],
                              pool_b: Sequence[str],
                              dedup: bool = True,
                              workers: int = 1,
                              chunksize: int = 64) -> list[dict]:
    """
    Full Cartesian A x B enumeration, optionally in parallel.

    Deduplication (by canonical product SMILES) is applied after the parallel
    combine step, preserving the first occurrence in (a, b) order.
    """
    n = resolve_workers(workers)
    pairs = [(i, a, j, b)
             for i, a in enumerate(pool_a)
             for j, b in enumerate(pool_b)]

    if n == 1 or len(pairs) <= 1:
        rows_iter = (_combine_pair(p) for p in pairs)
    else:
        with ProcessPoolExecutor(max_workers=n) as ex:
            rows_iter = list(ex.map(_combine_pair, pairs, chunksize=chunksize))

    rows, seen = [], set()
    for row in rows_iter:
        if row is None:
            continue
        if dedup:
            if row["product"] in seen:
                continue
            seen.add(row["product"])
        rows.append(row)
    return rows
