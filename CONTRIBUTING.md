# Contributing

Bug reports and pull requests are welcome. For bugs, include the Snakemake
version, the relevant rule log, the cluster job state, and a redacted
configuration file.

Before opening a pull request, run:

```bash
python -m unittest discover -s tests -v
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.test.yaml --dry-run --cores 1
```

Do not commit genome assemblies, sequence databases, generated results, Conda
environments, or credentials.
