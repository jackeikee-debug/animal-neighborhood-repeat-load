# Animal neighborhood repeat load

Code and reproducibility materials for repeat-load correspondence in homologous animal gene neighborhoods, including the three-hemichordate discovery panel and the 54-species orthology framework.

Author and correspondence information are provided in the manuscript. Licensing is pending author confirmation.
## Current scientific scope

The accepted manuscript reports recurring, heterogeneous correspondence with annotation, scale and membership dependence. The main expanded comparison retains34 finite estimates among45 candidates and all11 NA states. FullS35 and working intervals remain unavailable under their recorded gates. OG0000725 remains an unresolved154-gene/41-species exception.

## New synthetic calibration

Three prespecified mechanisms,100 outer replicates each,99 inner shifts. Rejection counts at.05/.10: S1 cyclic-group benchmark7/14, S2 nonstationary-variance stress6/11, S3 restricted non-group shifts2/10, each out of100. All six binomial95% intervals include nominal levels. None met the prespecified lower-bound inflation criterion. This finite diagnostic does not establish equivalence or calibrate real-data.001 tails. Full specification and S7 are included.

The numerical functions in physical_shift.py and physical_shift_root_numeric_gates.py are unchanged from the accepted implementation. The exported synthetic runner changes only output/import paths; every function AST matches the executed source. Historical stage guards, runtime binaries and real-data pipelines are separately documented dependencies. This compact repository is not a claim of standalone full numerical reproduction.

## Run a new independent synthetic namespace

Use a Python environment with NumPy and SciPy. Set PFLAVA_SYNTHETIC_ROOT to a new empty absolute output directory; create its replicate-claims, pilot and remainder subdirectories. From pilot run python /absolute/path/to/synthetic_zero.py pilot. From remainder run the same script with remainder. Each formal seed0..99 is exclusive per namespace; existing claim or result files stop a repeated invocation. The original project guard and its measured256MiB/1core budget are not automatically supplied by this standalone script.

## Evidence archive

Versioned v12 through v15 result archives and the new v16 diagnostic are being deposited in the restricted Zenodo draft https://zenodo.org/deposit/23723115. Upload receipts, accession and file checksums will be supplied after verification. Third-party Genome/GFF/RepeatMasker sources and full upstream assets remain identified external dependencies. No credentials are included. License selection is pending; no third-party data are relicensed by this repository.
