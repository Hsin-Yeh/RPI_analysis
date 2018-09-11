CXXFLAGS=-g -m64 -O2 -Wall
ROOTFLAGS=$(shell root-config --libs --cflags --glibs)
#SHAREFLAGS=-fPIC -c

All: makePlots	

#############################################################

#ped_check: main_ped.o hit_hits_class.so ped_check.o
#	g++ $^ -o ped_check $(CXXFLAGS) $(ROOTFLAGS)
#	mv *.o obj_file

#main_ped.o: main_ped.cc ped_check.h
#	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

#ped_check.o: ped_check.cxx ped_check.h hit_hits_class.h
#	g++ -c $(CXXFLAGS) $(ROOTFLAGS) ped_check.cxx -o $@

#############################################################

makePlots: main.o hit_hits_class.so makePlots.o
	g++ $^ -o makePlots $(CXXFLAGS) $(ROOTFLAGS)
	mv *.o obj_file

main.o: main_ped.cc makePlots.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

makePlots.o:makePlots.cxx makePlots.h hit_hits_class.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) makePlots.cxx -o $@

############################################################

clean:
	rm -f makePlots ped_check
	rm -f obj_file/*.o
