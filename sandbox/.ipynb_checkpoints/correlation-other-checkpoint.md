### For general:

https://xarray.pydata.org/en/stable/howdoi.html

### Causality analysis:

compute information flow myself first 
https://journals.aps.org/pre/pdf/10.1103/PhysRevE.90.052150

then try metrics
https://github.com/ZacKeskin/PyCausality
https://github.com/janekolszak/ifa


### For possibly other lines of approach:

tutorial
https://tutorial.xarray.dev/intermediate/01-high-level-computation-patterns.html

user guide
https://docs.xarray.dev/en/stable/user-guide/groupby.html

### For coarsen:

standardise dimension names
https://cf-xarray.readthedocs.io/en/latest/howtouse.html
https://xarray.pydata.org/en/stable/generated/xarray.Dataset.rename.html#xarray.Dataset.rename

correlation
https://docs.xarray.dev/en/stable/generated/xarray.corr.html

or align first then correlation
https://docs.xarray.dev/en/latest/generated/xarray.align.html

### For groupbybins:

stacked
https://docs.xarray.dev/en/stable/user-guide/groupby.html

groupbybins
https://docs.xarray.dev/en/stable/generated/xarray.DataArray.groupby_bins.html

or groupbybins then reindex
https://xarray.pydata.org/en/stable/generated/xarray.Dataset.reindex_like.html#xarray.Dataset.reindex_like

workarounds
https://discourse.pangeo.io/t/xarray-dataset-grouby-bins-without-squishing-other-dimensions/1692
https://github.com/pydata/xarray/issues/2488
https://stackoverflow.com/questions/40465026/groupby-bins-on-two-variables
https://github.com/pydata/xarray/issues/6610

### For xhistogram:

api
https://xhistogram.readthedocs.io/en/latest/api.html

tutorial
https://xhistogram.readthedocs.io/en/latest/tutorial.html

### For flox:

api
https://flox.readthedocs.io/en/latest/generated/flox.xarray.xarray_reduce.html

strategies
https://flox.readthedocs.io/en/latest/user-stories/climatology.html