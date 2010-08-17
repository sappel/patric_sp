MAKEFLAGS = $(MAKEOPTIONS)
MAKEOPTIONS = -j1  # use?; SP
CXX = mpiCC2
ifeq ($(MAKEFLAGS), debug)
  CXXOPTIMIZEFLAGS = -O0
  CXXFLAGS = -g3 $(CXXOPTIMIZEFLAGS)
else
  CXXOPTIMIZEFLAGS = -O3 -ffast-math
  CXXFLAGS = $(CXXOPTIMIZEFLAGS)
endif
INCLUDES = -I./include/ -I./tools/mpich/include -I./tools/obflib/include 
CCLINK = $(CC)
LD = $(CC)
LDOPTIONS = $(LDFLAGS)
LDFLAGS = -L. -L./tools/obflib/lib -lobflib -lmathtools -lstdc++   # -fexceptions 


bin/patric: obj/Main.o obj/SectorMap.o obj/Pic.o obj/TImpedance.o obj/Bump.o \
  ./tools/obflib/lib/libobflib.a
	$(CXX) -o $@ $^ $(LDOPTIONS) 

obj/Main.o: src/Main.cpp
	$(CXX) -c -o $@ $< $(INCLUDES) $(CXXFLAGS)

obj/SectorMap.o: src/SectorMap.cpp
	$(CXX) -c -o $@ $< $(INCLUDES) $(CXXFLAGS)

obj/Pic.o: src/Pic.cpp
	$(CXX) -c -o $@ $< $(INCLUDES) $(CXXFLAGS)

obj/TImpedance.o: src/TImpedance.cpp
	$(CXX) -c -o $@ $< $(INCLUDES) $(CXXFLAGS)

obj/Bump.o: src/Bump.cpp  # new: local orbit bump for injection; SP
	$(CXX) -c -o $@ $< $(INCLUDES) $(CXXFLAGS)

obj/Main.o: include/SectorMap.h	include/Pic.h

obj/SectorMap.o: include/SectorMap.h

obj/Pic.o: include/Pic.h

obj/TImpedance.o: include/TImpedance.h

obj/Bump.o: include/Bump.h  # SP

delfiles = obj/*  bin/patric python/odict.pyc *~ */*~ mad/madx01.eps mad/madx.ps
clean:
	@echo "Cleaning directory."
# delete files if exisiting and gives no error if not; SP
	@for i in $(delfiles); do [ -e "$$i" ] && rm "$$i" || 2>&1; done
