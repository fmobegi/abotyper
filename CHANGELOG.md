# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - [26-05-2026]

Update to 4.0.2 nf-core template.

### Changed

- aggregate_abo_reports.py: comprehensive A2 panel to prevent O1v alleles from triggering false A2 calls
- aggregate_abo_reports.py: lower c.1061delC / c.1032GA read threshold from 100 to 30 reads
- pysam_haploscan.py: new BAM/CIGAR-direct indel quantification module using pysam; resolves mpileup false-negative for rare variants. e.g IMM-26-23099 (32% del, was 0%)

## Fixed

- INDEL_DIAGNOSTIC pos687 corrected to pos685
- Validated against 51-sample cohort: all calls concordant with serology.

## v1.1.0 - [29-04-2026]

Update to 4.0.0 nf-core template.

## v1.1.0 - [02-04-2026]

Initial release of nf-core/abotyper, created with the [nf-core](https://nf-co.re/) template.

### Added

- Final ABO\*A1/A2 differentiation logic
- Testing conda docker and singularity profiles for requisite support after template sync

### Changed

- Bump version to 1.1.0 for first release
- Versions topics for local modules
- Updated inputs for samtools modules to match v.1.23.x

### Fixed

- Fixed MultiQC configuration errors to match topics
- Improved file staging for MultiQC compatibility for versions topics
- Updated `aggregate_abo_reports.py` to remove redundant variations
- Adrressed all pending review recommendation

## v1.0.0dev - [2025-09-01]

### Added

- Initial ABO\*A1/A2 differentiation algorithm (WORK IN PROGRESS)
- Metadata propagation for exonic information
- Improved meta.yaml documentation for local modules
- Updated documenting and staging of versions.yaml files for MultiQC
- Testing docker and singularity containers for requisite profiles support

### Changed

- Removed modules aliasing for cleaner code structure
- Improved staging for MultiQC files
- Updated MultiQC custom config to improve reporting
- Updated modules config file to match removed modules aliasing
- Moved index preparation to main workflow (streamlining processes)
- Upgraded all python scripts to class-style declaration
- Refactored all Python scripts to remove stale imports and improve maintainability

### Fixed

- Fixed MultiQC configuration errors
- Improved file staging for MultiQC compatibility
- Fixed issues with metadata propagation in the pipeline to reduce modules aliasing

### Known Issues

- MultiQC-1.30 treats the `\.` at the end of the `multiqc module command` as an illegal character causing failure. After long discussions and testing on different platforms, it was determined that this is a potential `AVX2` issue. Some python libraries may have been compiled with AVX2 instructions under containerised environments.

## [1.1.0] - [2025-07-31]

### New Features

- Initial release of ABO typing pipeline
- Comprehensive test suite with `test` profile support
- Debug mode compatibility for enhanced troubleshooting
- Updated usage documentation in `docs/usage.md`
- Updated output documentation in `docs/output.md`
- Tool citations and contributor acknowledgments in `README.md`

### Improvements

- Pipeline now fully compatible with nf-core standards
- Enhanced linting compliance (`nf-core pipelines lint`)
- Improved test dataset integration

### Bug Fixes

- All pipeline tests now pass successfully
- Debug mode warnings resolved
- Linting issues addressed

### Documentation Updates

- Updated `README.md` with new tool citations and contributors
- Comprehensive `docs/usage.md` documentation
- Detailed `docs/output.md` with pipeline outputs
- This `CHANGELOG.md` updated with all changes
