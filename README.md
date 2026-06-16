# Coordination Analysis

Computes per-atom coordination number and the radial distribution function (RDF).

## Install

```bash
vpm install @voltlabs/coordination-analysis
```

## CLI

```bash
coordination-analysis <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--cutoff <float>` | no | `3.2` | Cutoff radius for neighbor search. |
| `--rdf_bins <int>` | no | `500` | Number of bins for the RDF calculation. |
| `--threads <int>` | no | auto | Maximum worker threads. |
| `--help` | no | — | Print CLI help. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_coordination.parquet` | Coordination Analysis | — |
| `{output_base}_atoms.parquet` | Coordination Model | AtomisticExporter → glb |
| `{output_base}_rdf_chart.parquet` | RDF Chart | ChartExporter → chart-png |
| `{output_base}_rdf_histogram.parquet` | Partial RDF Histogram | — |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/coordination-analysis
