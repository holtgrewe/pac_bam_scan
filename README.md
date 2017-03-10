[![Build Status](https://travis-ci.org/holtgrewe/pac_bam_scan.svg?branch=master)](https://travis-ci.org/holtgrewe/pac_bam_scan)

# Ad-hoc PacBio BAM analysis

## Building

Look into `.travis-scripts/build.sh` for something working in the continuous integration.
Hopefully, the following is in sync:

```
# clone
$ git clone git@github.com:holtgrewe/pac_bam_scan.git
$ cd pac_bam_scan
$ git submodule init

# build
$ pushd htslib && sed -i.bak -e 's~^prefix.*~prefix=./inst~' Makefile && make install && popd
$ mkdir -p build && pushd build
$ export HTSLIB_ROOT=${PWD}/../htslib/inst
$ cmake -DCMAKE_CXX_COMPILER=$(which $CXX) -DCMAKE_C_COMPILER=$(which $CC) -DCMAKE_BUILD_TYPE=Release ..
$ make VERBOSE=1

# run
$ ./src/pac_bam_scan PATH_TO_FILE.bam >full_stats.txt
$ grep '^H\?READ' stats.txt | cut -f 3- >stats.txt
```

## Statistics

The Rmarkdown to build the statistics is in `reports/`
