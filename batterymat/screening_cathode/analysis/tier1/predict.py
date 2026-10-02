"""Inference harness for the archived tier-1 ALIGNN average-voltage model.

Loads ``best_model.pt`` + ``config.json`` from the archived Li_250 training
run and predicts the ALIGNN-FF average Li intercalation voltage (V) for
JARVIS ``atoms`` dicts, using exactly the graph construction the model was
trained with (k-nearest, cutoff 8.0 A, max_neighbors 12, CGCNN atom
features, canonized edges, line graph).

Usage (module):
    from predict import Tier1Predictor
    p = Tier1Predictor()
    v = p.predict_one(atoms_dict)            # float
    vs = p.predict_many(list_of_atoms_dicts) # list of float

Usage (CLI):
    python predict.py --sanity            # reproduce Test_results.json
    python predict.py --jids JVASP-2017 JVASP-42723   # from JARVIS cache
"""
import argparse
import json
import os
import pickle
import time

import numpy as np
import torch

ARCHIVE = "/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250"
MODEL_DIR = os.path.join(ARCHIVE, "voltage")
CACHE_PKL = os.path.join(ARCHIVE, "jarvis_dft3d_cache.pkl")


class Tier1Predictor:
    def __init__(self, model_dir=MODEL_DIR, checkpoint="best_model.pt", device="cpu"):
        from alignn.models.alignn_atomwise import ALIGNNAtomWise, ALIGNNAtomWiseConfig

        self.config = json.load(open(os.path.join(model_dir, "config.json")))
        torch.set_default_dtype(torch.float32)
        self.device = torch.device(device)
        mcfg = ALIGNNAtomWiseConfig(**self.config["model"])
        self.model = ALIGNNAtomWise(mcfg)
        state = torch.load(os.path.join(model_dir, checkpoint), map_location=self.device)
        if isinstance(state, dict) and "model" in state:
            state = state["model"]
        self.model.load_state_dict(state)
        self.model.to(self.device)
        self.model.eval()
        # graph settings straight from the training config
        self.cutoff = float(self.config["cutoff"])
        self.max_neighbors = int(self.config["max_neighbors"])
        self.use_canonize = bool(self.config["use_canonize"])
        self.cutoff_extra = float(self.config.get("cutoff_extra", 3.0))
        self.neighbor_strategy = self.config["neighbor_strategy"]
        self.atom_features = self.config["atom_features"]
        self.dtype = self.config.get("dtype", "float32")

    def graph(self, atoms):
        from jarvis.core.atoms import Atoms
        from alignn.graphs import Graph

        if isinstance(atoms, dict):
            atoms = Atoms.from_dict(atoms)
        g, lg = Graph.atom_dgl_multigraph(
            atoms,
            cutoff=self.cutoff,
            max_neighbors=self.max_neighbors,
            atom_features=self.atom_features,
            compute_line_graph=True,
            use_canonize=self.use_canonize,
            cutoff_extra=self.cutoff_extra,
            neighbor_strategy=self.neighbor_strategy,
            dtype=self.dtype,
        )
        lat = torch.tensor(atoms.lattice_mat).type(torch.get_default_dtype())
        return g, lg, lat

    @torch.no_grad()
    def predict_one(self, atoms):
        g, lg, lat = self.graph(atoms)
        out = self.model([g.to(self.device), lg.to(self.device), lat.to(self.device)])["out"]
        return float(out.detach().cpu().numpy().flatten()[0])

    def predict_many(self, atoms_list, verbose_every=200):
        preds = []
        t0 = time.time()
        for i, a in enumerate(atoms_list):
            try:
                preds.append(self.predict_one(a))
            except Exception as exc:  # keep going, mark failure
                print(f"  [warn] structure {i} failed: {exc}")
                preds.append(float("nan"))
            if verbose_every and (i + 1) % verbose_every == 0:
                el = time.time() - t0
                print(f"  {i+1}/{len(atoms_list)}  {el:.0f}s  ({el/(i+1):.3f} s/struct)", flush=True)
        return preds


def load_cache(path=CACHE_PKL):
    """Return {jid: entry} for the archived JARVIS dft_3d snapshot."""
    with open(path, "rb") as f:
        dft3d = pickle.load(f)
    return {e["jid"]: e for e in dft3d}


def sanity_check(n=None, tol=1e-3):
    """Reproduce archived Test_results.json predictions from id_prop.json."""
    test = json.load(open(os.path.join(MODEL_DIR, "Test_results.json")))
    idprop = {d["jid"]: d for d in json.load(open(os.path.join(ARCHIVE, "id_prop.json")))}
    p = Tier1Predictor()
    rows = test if n is None else test[:n]
    diffs = []
    t0 = time.time()
    for r in rows:
        v = p.predict_one(idprop[r["id"]]["atoms"])
        diffs.append(abs(v - r["pred_out"]))
        print(f"{r['id']:14s} archived={r['pred_out']:+.5f} recomputed={v:+.5f} |d|={diffs[-1]:.2e} label={r['target_out'][0]:+.3f}")
    diffs = np.array(diffs)
    print(f"\n{len(rows)} test structures in {time.time()-t0:.1f}s; max|d|={diffs.max():.2e}, mean|d|={diffs.mean():.2e}, n>{tol}: {(diffs>tol).sum()}")
    return diffs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sanity", action="store_true")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--jids", nargs="*")
    a = ap.parse_args()
    if a.sanity:
        sanity_check(a.n)
    if a.jids:
        cache = load_cache()
        p = Tier1Predictor()
        for j in a.jids:
            print(j, cache[j]["formula"], f"{p.predict_one(cache[j]['atoms']):.4f} V")
