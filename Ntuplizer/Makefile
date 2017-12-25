CXXFLAGS=-g -m64 -O2 -Wall
ROOTFLAGS=$(shell root-config --libs --cflags --glibs) -lMathMore
SHAREFLAGS=-fPIC -c

Ntuplizer:  CHMapping.o  main.o  hit_hits_class.so
	g++ $^ -o Ntuplizer $(CXXFLAGS) $(ROOTFLAGS)
	cp *.so *.pcm ../Analysis
	mv *.o  *.so  *.cxx obj_file

CHMapping.o: CHMapping.cc CHMapping.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) CHMapping.cc -o  CHMapping.o

main.o: decode_main.cc CHMapping.h hit_hits_class.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) decode_main.cc -o  main.o

hit_hits_class.cxx: hit_hits_class.cc hit_hits_class.h linkdef.h
	rootcint -f  $@ -c hit_hits_class.cc hit_hits_class.h linkdef.h

hit_hits_class.so: hit_hits_class.cxx
	g++ $(CXXFLAGS) $(ROOTFLAGS) $(SHAREFLAGS) -o  $@ hit_hits_class.cxx

clean:
	rm  obj_file/*.o  obj_file/*.so  obj_file/*.cxx *.pcm
	rm Ntuplizer
