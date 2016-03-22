SAMTOOLS_DIR=samtools
LIBBAM=bam
osname := $(shell uname -s 2>/dev/null || echo "NA")
ifeq ($(osname),Darwin)
    LIBBAM=bammac
endif

all: rnapileup filter_pileup rnapileup2mismatchbed 

#samtoolssrc:
#	$(MAKE) -C samtools

rnapileup: rnapileup.cpp
	g++ -O2 -L $(SAMTOOLS_DIR) -I $(SAMTOOLS_DIR) $? -o $@ -l$(LIBBAM) -lz

filter_pileup: filter_pileup.cpp
	g++ -O2 -L $(SAMTOOLS_DIR) -I $(SAMTOOLS_DIR) $? -o $@ -l$(LIBBAM) -lz

rnapileup2mismatchbed: rnapileup2mismatchbed.cpp
	g++ -O2 $? -o $@

clean:
	rm -f rnapileup filter_pileup rnapileup2mismatchbed


