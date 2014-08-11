# Requirements

* Operating system: Linux. This packeage relies on some linux utilities for text processing (sed, cut etc).
It's possible to install these tools on Windows through MinGW or Cygwin, but I haven't tried it.
* R packages:
    * devtools
    * ggplot2

# Installation

    require(devtools)
    Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
    install_github("bedcollr", username="kindlychung")

# Changes

* added R.utils dependency
* contrastplot now working under new design
* readout now working under the new design
* xxx_shift_0000.bed now linked to original xxx.bed
