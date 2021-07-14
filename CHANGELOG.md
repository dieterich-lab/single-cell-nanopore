# Change Log
Notable changes to scNapBar will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Changed
- Updated **NanoSim** version 3.0.1 with fix (`environment.yaml`, documentation).
- Modify sort in `align_longreads`, and pipe output to BAM directly.
- Minor changes: move hard coded fields to config file, adjust `cluster.json`, *etc.*
- Updated documentation.

### Added
- File ./gitignore to ignore build/temporary files.
- Seed to `sim_reads` (NanoSim `simulator.py`), and `shuf` in `build_genome`.

### Fixed
- Delete temporary files.
- Use threads.

### Removed
- Rule `unzip_fq`.

## [1.0.0] Undocumented
