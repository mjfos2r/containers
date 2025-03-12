MODULES = align-tools chopper canu basic hgtector kraken2 mummer4 nanoplot samtools

all:
	for dir in $(MODULES); do \
		echo begin $$dir; \
		(cd $$dir; ${MAKE}); \
		echo end $$dir; \
		echo ==========; \
	done