# CoordinationAnalysis

`CoordinationAnalysis` computes coordination statistics and radial distribution data.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- CoordinationAnalysis
```

## CLI

Usage:

```bash
coordination-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--cutoff <float>` | No | Cutoff radius for neighbor search. | `3.2` |
| `--rdfBins <int>` | No | Number of bins for RDF calculation. | `500` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |
