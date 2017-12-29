CXXFLAGS=-g -std=c++11 -m64 -O2 -Wall
ROOTFLAGS=$(shell root-config --libs --cflags --glibs)
SHAREFLAGS=-fPIC -c

all: Ntuplizer COPY

COPY:
	$(eval pcm_exist=$(wildcard *.pcm))
ifneq ($(strip $(pcm_exist)), )  #ROOT6 has different output from rootcint(.pcm)
	cp *.so *.pcm ../Analysis
	mv *.o  *.so  *.cxx obj_file
else                             #Bofore ROOT6 version
	cp *.so ../Analysis
	mv *.o  *.so *.cxx obj_file
endif

Ntuplizer:  CHMapping.o  main.o  hit_hits_class.so
	g++ $^ -o Ntuplizer $(CXXFLAGS) $(ROOTFLAGS)

CHMapping.o: CHMapping.cc CHMapping.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) CHMapping.cc -o  CHMapping.o

main.o: decode_main.cc CHMapping.h hit_hits_class.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) decode_main.cc -o  main.o

hit_hits_classdict.cxx: hit_hits_class.cc hit_hits_class.h linkdef.h
	rootcint -f  $@ -c hit_hits_class.cc hit_hits_class.h linkdef.h

hit_hits_class.so: hit_hits_classdict.cxx
	g++ $(CXXFLAGS) $(ROOTFLAGS) $(SHAREFLAGS) -o  $@ hit_hits_classdict.cxx

clean:
	rm -f obj_file/*.o  obj_file/*.so  obj_file/*.cxx *dict*
	rm -f Ntuplizer
