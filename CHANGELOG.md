# Change Log
Notable changes to scNapBar will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Fixed
- Modify seed for `shuf` in `build_genome`. 

## [1.1.0] 15.10.2021

### Changed
- Update **NanoSim** to latest version (environment.yaml), simplify install.
- Modify sort in `align_longreads`, and pipe output to BAM directly.
- Minor changes: move hard coded fields to config file, adjust `cluster.json`, *etc.*
- Updated documentation.
- Upgraded minimap2 to 2.21-r1071.

### Added
- File ./gitignore to ignore build/temporary files.
- Seed to `sim_reads` (NanoSim `simulator.py`), and `shuf` in `build_genome`.
- I/O option to minimap2 call in `align_longreads`.

### Fixed
- Delete temporary files.
- Use threads. TODO: see [Control the number of cores/threads per rule](https://github.com/dieterich-lab/single-cell-nanopore/issues/14).
- Temporary fix in `sim_reads` due to [Breaking changes in format output from 3.0.0-beta](https://github.com/bcgsc/NanoSim/issues/132). TODO: modify pipeline.

### Removed
- Rule `unzip_fq`.

## [1.0.0] Undocumented
