# ebdetect
[![github](https://img.shields.io/badge/GitHub-grviss%2Febdetect-blue.svg?style=flat)](https://github.com/grviss/ebdetect)
[![ADS](https://img.shields.io/badge/ADS-arXiv190107975V-red.svg)](https://ui.adsabs.harvard.edu/abs/2019arXiv190107975V/abstract)
[![arxiv](https://img.shields.io/badge/arxiv-1901.07975-orange.svg?style=flat)](https://arxiv.org/abs/1901.07975)

`EBDETECT` is a solar feature detection and tracking code developed to identify Ellerman and UV bursts based on intensity thresholding, as well as size, continuity and lifetime constraints. The procedure and results from detection in SST H&alpha; and SDO/AIA 1700Å data are described in [Vissers et al. (2019)](https://arxiv.org/abs/1901.07975) (see above ADS and arXiv links).

The code is provided "as-is". Bug reports are appreciated and encouraged, but dedicated support will only be provided within a scientific collaboration.

## Running EBDETECT
The minimal input required for `EBDETECT` to run is a plain-text configuration file, wherein the data file (currently only legacy SST/La Palma &mdash; i.e. CRISPEX-ready &mdash; format allowed) that detection will be performed on is specified.
`EBDETECT` parses the input filename(s), detection settings and I/O options from the configuration file (see `ebdetect_config.txt` for a template). 

The detection process consists of four steps:
1. Intensity thresholding per frame. Threshold defaults to mean intensity, but can be adjusted using either `SDEV_MULT_CONSTRAINT`, `MEAN_MULT_CONSTRAINT` or `INTENSITY_CONSTRAINT`. 
2. Size constraint (lower or lower and upper limits). Defaults to `[0,nx*ny]` range, but can be adjusted using `SIZE_CONSTRAINT` to a scalar or 2-element array.
3. Overlap and continuity (minimum spatial overlap between frames and maximum skipped frames). Defaults to minimum overlap of 1 pixel (adjusted using `OVERLAP_CONSTRAINT`) and 0 skipped frames (adjusted using `T_SKIP_CONSTRAINT`).
4. Lifetime constraint (lower or lower and upper limits). Defaults to `[0,nt]` range, but can be adjusted using `LIFETIME_CONSTRAINT` to a scalar or 2-element array.

## Configuration file
### Required keywords
* `INPUTFILE`: full filename (without path). Currently restricted to being a legacy SST/La Palma format cube (i.e. data order either `[nx, ny, nt]` or `[nx, ny, nw, nt]`, with `nw` and `nt` the size of the wavelength and time dimensions, respectively).
* `NW`: size of wavelength dimension in `INPUTFILE`. Defaults to 1.

### Optional keywords
#### Intensity thresholding
* `SDEV_MULT_CONSTRAINT`: multiplier factor to determine the threshold as fraction of sigma above mean. Threshold gets calculated as: `threshold = I_mean + SDEV_MULT_CONSTRAINT * I_stdev`
* `MEAN_MULT_CONSTRAINT`: multiplier factor to determine threshold as fraction of mean. Overrides setting of `SDEV_MULT_CONSTRAINT`. Threshold gets calculated as: `threshold = MEAN_MULT_CONSTRAINT * I_mean`
* `INTENSITY_CONSTRAINT`: intensity threshold in counts. Takes precedence over both `SDEV_MULT_CONSTRAINT` and `MEAN_MULT_CONSTRAINT`.  Threshold gets calculated as: `threshold = INTENSITY_CONSTRAINT`

## Output
`EBDETECT` outputs an IDL save file with a pointer structure containing the detections. When `WRITE_MASK=1` is set in the configuration file (default behaviour), a data cube with boolean masks highlighting the detected pixels for every frame will also be saved.

Output can be saved at three stages in the detection process:
* After the intensity thresholding (detection step #1; `WRITE_DETECT_INIT=1`)
* After applying size constraint(s), and checking for detection overlap and merging/splitting (detection step #3; `WRITE_DETECT_OVERLAP=1`)
* After detection grouping and applying lifetime constraint(s), i.e. the final detection results (detection step #4; `WRITE_DETECT_FINAL=1`)

By default these keywords are set and output will be saved after all three of the above stages. Note that `WRITE_MASK` combines with the above keywords and cannot be set independently for these save points (e.g. if `WRITE_DETECT_OVERLAP=0` and `WRITE_MASK=1` then `EBDETECT` will output both an IDL save file and a boolean mask file after intensity thresholding and after compiling the final results, but save no output after the overlap and continuity check).
