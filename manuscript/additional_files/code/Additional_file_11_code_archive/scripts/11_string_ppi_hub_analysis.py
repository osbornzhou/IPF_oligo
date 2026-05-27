#!/usr/bin/env python
"""
Fetch STRING PPI networks and compute hub genes with triple QC.

Primary network:
  - robust_mrna_strict: robust IPF mRNA candidates from discovery-validation analysis

Secondary network:
  - mirna_mrna_axis_targets_all: target genes from robust negative miRNA-mRNA axes

Outputs:
  - results/ppi_network/string_ppi_edges_*.csv
  - results/ppi_network/string_ppi_nodes_*.csv
  - results/ppi_network/string_ppi_hub_genes_*.csv
  - results/ppi_network/string_ppi_triple_qc.csv
  - results/ppi_network/string_ppi_summary.xlsx
"""

from __future__ import annotations

import csv
import io
import math
import time
import urllib.parse
import urllib.request
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
ROBUST_DIR = PROJECT_DIR / "results" / "robust_candidates"
AXES_DIR = PROJECT_DIR / "results" / "mirna_mrna_axes"
OUTPUT_DIR = PROJECT_DIR / "results" / "ppi_network"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

STRING_API = "https://string-db.org/api"
SPECIES = 9606
CALLER_IDENTITY = "ipf_oligo_ml_bmc_genomics_project"
REQUIRED_SCORES = {
    "medium_confidence": 400,
    "high_confidence": 700,
}


def normalize_gene(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.upper()
    return text if text and text != "NA" else ""


def string_api_post(method: str, params: dict[str, object], retries: int = 3) -> str:
    url = f"{STRING_API}/tsv/{method}"
    data = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(url, data=data, method="POST")
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    last_error = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"STRING API request failed for {method}: {last_error}")


def string_api_get(method: str, params: dict[str, object] | None = None) -> str:
    params = params or {}
    url = f"{STRING_API}/tsv/{method}"
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def parse_tsv(text: str) -> pd.DataFrame:
    if not text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(text), sep="\t")


def load_gene_sets() -> dict[str, list[str]]:
    robust = pd.read_csv(ROBUST_DIR / "robust_mrna_candidates_strict.csv")
    axes = pd.read_csv(AXES_DIR / "robust_mirna_mrna_negative_axes_mirtarbase.csv")

    robust_genes = sorted({normalize_gene(x) for x in robust["standard_feature_id"] if normalize_gene(x)})
    axis_targets = sorted({normalize_gene(x) for x in axes["target_gene"] if normalize_gene(x)})

    return {
        "robust_mrna_strict": robust_genes,
        "mirna_mrna_axis_targets_all": axis_targets,
    }


def get_string_version() -> str:
    try:
        text = string_api_get("version")
        return " ".join(text.strip().split())
    except Exception as exc:  # noqa: BLE001
        return f"STRING version query failed: {exc}"


def map_genes(gene_set_name: str, genes: list[str]) -> pd.DataFrame:
    text = string_api_post(
        "get_string_ids",
        {
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "limit": 1,
            "echo_query": 1,
            "caller_identity": CALLER_IDENTITY,
        },
    )
    mapped = parse_tsv(text)
    if mapped.empty:
        return pd.DataFrame()
    mapped["gene_set"] = gene_set_name
    mapped["queryItem"] = mapped["queryItem"].map(normalize_gene)
    mapped["preferredName"] = mapped["preferredName"].map(normalize_gene)
    mapped.to_csv(OUTPUT_DIR / f"string_mapping_{gene_set_name}.csv", index=False, encoding="utf-8-sig")
    return mapped


def fetch_network(gene_set_name: str, genes: list[str], score_label: str, required_score: int) -> pd.DataFrame:
    text = string_api_post(
        "network",
        {
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "required_score": required_score,
            "network_type": "functional",
            "add_nodes": 0,
            "caller_identity": CALLER_IDENTITY,
        },
    )
    edges = parse_tsv(text)
    if edges.empty:
        return pd.DataFrame()

    edges["gene_set"] = gene_set_name
    edges["score_label"] = score_label
    for col in ["preferredName_A", "preferredName_B"]:
        edges[col] = edges[col].map(normalize_gene)
    if "score" in edges.columns:
        edges["combined_score"] = pd.to_numeric(edges["score"], errors="coerce")
    elif "combined_score" in edges.columns:
        edges["combined_score"] = pd.to_numeric(edges["combined_score"], errors="coerce")
    else:
        edges["combined_score"] = math.nan
    edges = edges[edges["preferredName_A"] != edges["preferredName_B"]].copy()
    edges = edges.sort_values("combined_score", ascending=False)
    edges.to_csv(OUTPUT_DIR / f"string_ppi_edges_{gene_set_name}_{score_label}.csv", index=False, encoding="utf-8-sig")
    return edges


def build_graph(edges: pd.DataFrame) -> dict[str, dict[str, float]]:
    graph: dict[str, dict[str, float]] = defaultdict(dict)
    for _, row in edges.iterrows():
        a = row["preferredName_A"]
        b = row["preferredName_B"]
        score = float(row["combined_score"])
        if not a or not b or a == b or not math.isfinite(score):
            continue
        graph[a][b] = max(score, graph[a].get(b, 0.0))
        graph[b][a] = max(score, graph[b].get(a, 0.0))
    return graph


def connected_components(graph: dict[str, dict[str, float]]) -> list[set[str]]:
    seen = set()
    components = []
    for node in graph:
        if node in seen:
            continue
        comp = set()
        queue = deque([node])
        seen.add(node)
        while queue:
            current = queue.popleft()
            comp.add(current)
            for neighbor in graph[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(comp)
    return components


def closeness_centrality(graph: dict[str, dict[str, float]]) -> dict[str, float]:
    centrality = {}
    nodes = list(graph)
    n = len(nodes)
    for source in nodes:
        distances = {source: 0}
        queue = deque([source])
        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor not in distances:
                    distances[neighbor] = distances[current] + 1
                    queue.append(neighbor)
        if len(distances) <= 1:
            centrality[source] = 0.0
        else:
            centrality[source] = (len(distances) - 1) / sum(distances.values())
            if len(distances) < n:
                centrality[source] *= (len(distances) - 1) / (n - 1)
    return centrality


def betweenness_centrality(graph: dict[str, dict[str, float]]) -> dict[str, float]:
    """Unweighted Brandes betweenness centrality."""
    nodes = list(graph)
    cb = {v: 0.0 for v in nodes}
    for s in nodes:
        stack = []
        pred = {w: [] for w in nodes}
        sigma = dict.fromkeys(nodes, 0.0)
        dist = dict.fromkeys(nodes, -1)
        sigma[s] = 1.0
        dist[s] = 0
        queue = deque([s])
        while queue:
            v = queue.popleft()
            stack.append(v)
            for w in graph[v]:
                if dist[w] < 0:
                    queue.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = dict.fromkeys(nodes, 0.0)
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w] != 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                cb[w] += delta[w]
    if len(nodes) > 2:
        scale = 1.0 / ((len(nodes) - 1) * (len(nodes) - 2))
        cb = {node: value * scale for node, value in cb.items()}
    return cb


def clustering_coefficient(graph: dict[str, dict[str, float]]) -> dict[str, float]:
    coeff = {}
    for node, neighbors in graph.items():
        neigh = list(neighbors)
        k = len(neigh)
        if k < 2:
            coeff[node] = 0.0
            continue
        links = 0
        for i in range(k):
            for j in range(i + 1, k):
                if neigh[j] in graph[neigh[i]]:
                    links += 1
        coeff[node] = (2 * links) / (k * (k - 1))
    return coeff


def compute_node_metrics(gene_set_name: str, score_label: str, input_genes: list[str], mapped: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    graph = build_graph(edges)
    components = connected_components(graph)
    component_id = {}
    component_size = {}
    for idx, comp in enumerate(sorted(components, key=len, reverse=True), start=1):
        for node in comp:
            component_id[node] = idx
            component_size[node] = len(comp)

    closeness = closeness_centrality(graph)
    betweenness = betweenness_centrality(graph)
    clustering = clustering_coefficient(graph)

    mapped_lookup = mapped.drop_duplicates("preferredName").set_index("preferredName").to_dict(orient="index") if not mapped.empty else {}
    rows = []
    all_nodes = sorted(set(input_genes) | set(graph))
    for node in all_nodes:
        neighbors = graph.get(node, {})
        weighted_degree = sum(neighbors.values())
        degree = len(neighbors)
        rows.append(
            {
                "gene_set": gene_set_name,
                "score_label": score_label,
                "gene_symbol": node,
                "in_input": node in set(input_genes),
                "mapped_to_string": node in mapped_lookup,
                "string_id": mapped_lookup.get(node, {}).get("stringId", ""),
                "degree": degree,
                "weighted_degree": weighted_degree,
                "betweenness": betweenness.get(node, 0.0),
                "closeness": closeness.get(node, 0.0),
                "clustering_coefficient": clustering.get(node, 0.0),
                "component_id": component_id.get(node, 0),
                "component_size": component_size.get(node, 1 if node in set(input_genes) else 0),
            }
        )
    nodes = pd.DataFrame(rows)
    for metric in ["degree", "weighted_degree", "betweenness", "closeness"]:
        max_value = nodes[metric].max()
        nodes[f"{metric}_norm"] = nodes[metric] / max_value if max_value and max_value > 0 else 0
    nodes["hub_score"] = (
        nodes["degree_norm"]
        + nodes["weighted_degree_norm"]
        + nodes["betweenness_norm"]
        + nodes["closeness_norm"]
    )
    nodes = nodes.sort_values(["hub_score", "degree", "weighted_degree"], ascending=False)
    nodes["hub_rank"] = range(1, len(nodes) + 1)
    nodes["hub_top20"] = nodes["hub_rank"] <= 20
    nodes.to_csv(OUTPUT_DIR / f"string_ppi_nodes_{gene_set_name}_{score_label}.csv", index=False, encoding="utf-8-sig")
    nodes[nodes["hub_top20"]].to_csv(OUTPUT_DIR / f"string_ppi_hub_genes_{gene_set_name}_{score_label}.csv", index=False, encoding="utf-8-sig")
    return nodes


def main() -> None:
    gene_sets = load_gene_sets()
    string_version = get_string_version()
    qc_rows = []
    workbook_sheets = {}

    for gene_set_name, genes in gene_sets.items():
        mapped = map_genes(gene_set_name, genes)
        mapped_genes = sorted(set(mapped["preferredName"]) & set(genes)) if not mapped.empty else []
        mapped_preferred_genes = sorted(set(mapped["preferredName"])) if not mapped.empty else []
        mapping_rate = len(mapped_genes) / max(len(genes), 1)

        for score_label, required_score in REQUIRED_SCORES.items():
            edges = fetch_network(gene_set_name, genes, score_label, required_score)
            nodes = compute_node_metrics(gene_set_name, score_label, mapped_preferred_genes, mapped, edges)

            node_genes = set(nodes.loc[nodes["degree"] > 0, "gene_symbol"])
            input_gene_set = set(genes)
            edge_nodes = set(edges["preferredName_A"]) | set(edges["preferredName_B"]) if not edges.empty else set()
            score_ok = edges["combined_score"].between(0, 1).all() if not edges.empty else False
            no_self_loops = (edges["preferredName_A"] != edges["preferredName_B"]).all() if not edges.empty else False
            edges_subset_input = edge_nodes.issubset(set(mapped_preferred_genes))

            qc_rows.append(
                {
                    "gene_set": gene_set_name,
                    "score_label": score_label,
                    "required_score": required_score,
                    "input_genes": len(genes),
                    "mapped_genes": len(mapped_genes),
                    "mapping_rate": round(mapping_rate, 6),
                    "edges": len(edges),
                    "network_nodes_with_edges": len(node_genes),
                    "connected_components": len(connected_components(build_graph(edges))) if not edges.empty else 0,
                    "largest_component_size": int(nodes["component_size"].max()) if not nodes.empty else 0,
                    "qc1_input_pass": len(genes) >= 5 and len(set(genes)) == len(genes),
                    "qc2_string_mapping_pass": mapping_rate >= 0.80 and len(mapped_genes) >= 5,
                    "qc3_network_integrity_pass": len(edges) > 0 and score_ok and no_self_loops and edges_subset_input,
                    "triple_qc_pass": (len(genes) >= 5 and len(set(genes)) == len(genes))
                    and (mapping_rate >= 0.80 and len(mapped_genes) >= 5)
                    and (len(edges) > 0 and score_ok and no_self_loops and edges_subset_input),
                    "string_version": string_version,
                    "organism": "Homo sapiens",
                    "species_taxon_id": SPECIES,
                }
            )

            workbook_sheets[f"{gene_set_name[:18]}_{score_label[:4]}_hub"] = nodes.head(50)
            workbook_sheets[f"{gene_set_name[:18]}_{score_label[:4]}_edges"] = edges.head(5000)

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUTPUT_DIR / "string_ppi_triple_qc.csv", index=False, encoding="utf-8-sig")
    workbook_sheets = {"triple_qc": qc, **workbook_sheets}

    with pd.ExcelWriter(OUTPUT_DIR / "string_ppi_summary.xlsx", engine="openpyxl") as writer:
        for sheet_name, df in workbook_sheets.items():
            safe_sheet = sheet_name[:31]
            df.to_excel(writer, sheet_name=safe_sheet, index=False)
            writer.book[safe_sheet].freeze_panes = "A2"

    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
