target=simulation
OO=mathfunc.o faddeeva.o acc_model.o beam.o beambeam.o ${target}.o beamu.o
${target}:${OO}
	mpic++ -o ${target} ${OO}
${target}.o:../beambeam.h ../beam.h ../acc_model.h ../faddeeva.h ../mathfunc.h ../${target}.cpp
	mpic++ -c -O3 ../${target}.cpp -std=c++11
mathfunc.o:../mathfunc.h ../mathfunc.cpp
	mpic++ -c -O3 ../mathfunc.cpp -std=c++11
faddeeva.o:../faddeeva.h ../faddeeva.cpp
	mpic++ -c -O3 ../faddeeva.cpp -std=c++11
acc_model.o:../acc_model.h ../mathfunc.h ../acc_model.cpp
	mpic++ -c -O3 ../acc_model.cpp -std=c++11
beam.o:../beam.h ../acc_model.h ../faddeeva.h ../mathfunc.h ../beam.cpp
	mpic++ -c -O3 ../beam.cpp -std=c++11
beambeam.o:../beambeam.h ../beam.h ../acc_model.h ../faddeeva.h ../mathfunc.h ../beambeam.cpp
	mpic++ -c -O3 ../beambeam.cpp -std=c++11
beamu.o:../beam.h ../beamu.cpp
	mpic++ -c -O3 ../beamu.cpp -std=c++11
.PHONY: clean
clean:
	rm -f *.o main
