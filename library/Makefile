DLibrary.so:	DData.cxx DParticle.cxx DEvent.cxx dict.cxx
	g++ -shared -o $@ `root-config --cflags --libs` -fPIC $^

dict.cxx:	DData.h DParticle.h DEvent.h LinkDef.h
	rootcling -f $@ -rml DLibrary.so -p $^