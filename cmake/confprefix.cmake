set(ST_LATEX_DIR "/usr/bin/"
	CACHE STRING
	"LaTeX default prefix (/usr/bin/)"
	)
set(ST_METIS_DIR $ENV{METIS}
	CACHE STRING
	"METIS should be compiled with same compiler GNU|Intel."
	)
message(STATUS "METIS locates ${ST_METIS_DIR}")
