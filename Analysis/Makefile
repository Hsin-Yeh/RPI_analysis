CXXFLAGS=-g -m64 -O2 -Wall
ROOTFLAGS=$(shell root-config --libs --cflags --glibs)
#SHAREFLAGS=-fPIC -c

All: makePlots	

makePlots: main.o hit_hits_class.so makePlots.o
	g++ $^ -o makePlots $(CXXFLAGS) $(ROOTFLAGS)
	mv *.o obj_file

main.o: main.cc makePlots.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

makePlots.o:makePlots.cxx makePlots.h hit_hits_class.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) makePlots.cxx -o $@

############################################################

clean:
	rm -f makePlots
	rm -f obj_file/*.o
