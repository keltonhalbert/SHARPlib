qc
==

QC builds, the default, treat a value that is ``sharp::MISSING`` or NaN as
missing, and every routine that skips or bridges missing data checks it with
``sharp::is_missing``. Builds configured with ``-DNO_QC=ON`` skip these checks.

.. doxygenfunction:: sharp::is_missing
