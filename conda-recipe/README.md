# conda-forge recipe

This directory mirrors the recipe submitted to
[conda-forge/staged-recipes](https://github.com/conda-forge/staged-recipes) as
`recipes/ncountr/recipe.yaml`. It is kept here so the recipe is versioned
alongside the package it builds.

It uses the [CEP 13][cep13] v1 recipe format (`recipe.yaml`), which is what
staged-recipes currently expects; the older conda-build `meta.yaml` format is
deprecated there.

When cutting a release, update `context.version` and `source.sha256` together.
The checksum is the one PyPI publishes for the sdist:

```bash
curl -sL "https://pypi.org/pypi/ncountr/<version>/json" \
  | python -c "import json,sys; print(next(u['digests']['sha256'] for u in json.load(sys.stdin)['urls'] if u['packagetype']=='sdist'))"
```

[cep13]: https://github.com/conda/ceps/blob/main/cep-13.md
