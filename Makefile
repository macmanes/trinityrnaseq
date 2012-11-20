###################################################################
#
# The default compiler is GNU gcc/g++.
# Run
#  make TRINITY_COMPILER=intel
# to build Inchworm and Chrysalis with the Intel compiler.
#

ifeq ($(TRINITY_COMPILER),intel)
INCHWORM_CONFIGURE_FLAGS = CXX=icpc CXXFLAGS="-fast"
CHRYSALIS_MAKE_FLAGS = COMPILER=icpc
else
override TRINITY_COMPILER=gnu
endif


all: 
	@echo Using $(TRINITY_COMPILER) compiler for Inchworm and Chrysalis
	cd Inchworm && (test -e configure || autoreconf) \
                && ./configure --prefix=`pwd` $(INCHWORM_CONFIGURE_FLAGS) && $(MAKE) install
	cd Chrysalis && $(MAKE) UNSUPPORTED=yes $(CHRYSALIS_MAKE_FLAGS)
	cd trinity-plugins/rsem && $(MAKE) 
	cd trinity-plugins/jellyfish && ./configure --prefix=`pwd` && $(MAKE)
	cd trinity-plugins/fastool && $(MAKE)
	cd trinity-plugins/slclust && $(MAKE)
	cd trinity-plugins/coreutils && ./build_parallel_sort.sh
	cd trinity-plugins/collectl && ./build_collectl.sh

clean:
	cd Inchworm && make clean
	cd Chrysalis && $(MAKE) clean UNSUPPORTED=yes
	cd trinity-plugins/rsem && $(MAKE) clean
	cd trinity-plugins/jellyfish && $(MAKE) clean 
	cd trinity-plugins/fastool && $(MAKE) clean
	cd trinity-plugins/slclust && $(MAKE) clean
	cd trinity-plugins/collectl && rm -rf bin doc man
	cd trinity-plugins/coreutils && rm -rf bin


###################################################################


