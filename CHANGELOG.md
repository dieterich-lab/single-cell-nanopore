# Change Log
Notable changes to scNapBar will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Changed
- Updated **NanoSim** latest master branch with fix (`environment.yaml`, documentation). TODO: fix release as soon as available.
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
- Use threads.
- Temporary fix in `sim_reads` due to https://github.com/bcgsc/NanoSim/issues/132. TODO: find solution.

### Removed
- Rule `unzip_fq`.

## [1.0.0] Undocumented
