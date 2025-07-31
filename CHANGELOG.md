# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/abotyper, created with the [nf-core](https://nf-co.re/) template.

### Added

- Initial ABO*A1/A2 differentiation algorithm (WORK IN PROGRESS)
- Metadata propagation for exonic information
- Improved meta.yaml documentation for local modules
- Updated documenting and staging of versions.yaml files for MultiQC

### Changed

- Removed modules aliasing for cleaner code structure
- Improved staging for MultiQC files
- Updated MultiQC custom config to improve reporting
- Updated modules config file to match removed modules aliasing
- Moved index preparation to main workflow (streamlining processes)

### Fixed

- Fixed MultiQC configuration syntax errors
- Improved file staging for MultiQC compatibility

### Known Issues

- TODO: Explore why MultiQC-1.30 treats the `\.` at the end of command as an illegal character causing failure

## [1.0.0] - 2025-07-31

### Added

- Initial release of ABO typing pipeline
