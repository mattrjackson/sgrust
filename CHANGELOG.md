# Change log
<!-- Template:
## [version] — YYYY-MM-DD

### Security
### Added
### Changed
### Deprecated
### Removed
### Fixed
-->

## [0.2.1] — 2024-12-20

### Security
### Added
Created a new `refine_parallel` method for the `SparseGrid` trait. The results in values being filled for the refined nodes in parallel.
### Changed
`RefinementFunctor` must now implement Send + Sync.
### Deprecated
### Removed
Removed the useless sorting in the `refine` function. Presently we do not limit the number of nodes to refine in a given iteration.
### Fixed

## [0.2.0] — 2024-12-16

### Security
### Added
### Changed
Improved performance of `CombinationGrid` interpolation by ~40%. It is now on par with the performance of TASMANIAN.
### Deprecated
### Removed
### Fixed

## [0.1.1] — 2024-12-14

### Security
### Added
### Changed
### Deprecated
### Removed
Removed unused tests file.
### Fixed

## [0.1.0] — 2024-12-14

### Security
### Added
Initial commit of base capabilities
### Changed
### Deprecated
### Removed
### Fixed