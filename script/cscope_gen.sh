#!/bin/bash
find . -name 'st*.fi' \
	-o -name 'st*.for'\
	-o -name 'st*.[fF]'\
	-o -name 'st*.FOR' \
	-o -name 'st*.f90' \
> cscope.files

# -b : just build
# -q : create inverted index
cscope -b -q
