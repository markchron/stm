#!/bin/bash
find . -name '*.fi' \
	-o -name '*.for'\
	-o -name '*.[fF]'\
	-o -name '*.FOR' \
	-o -name '*.f90' \
> cscope.files

# -b : just build
# -q : create inverted index
cscope -b -q
