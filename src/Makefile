CXXFLAGS=-O3 -std=c++11
PREFIX=/usr/local/
HEADERS=$(PREFIX)include/
BINARIES=$(PREFIX)bin/
MANPAGES=$(PREFIX)man/man1/
prefix=
rasa=-lfreesasa -DFREESASA
build=


.PHONY: default fragments all install clean clean-all

default: inspire

fragments: $(build)$(prefix)assign $(build)$(prefix)classify $(build)$(prefix)exclude $(build)$(prefix)features $(build)$(prefix)filter $(build)$(prefix)fingerprints $(build)$(prefix)index $(build)$(prefix)merge $(build)$(prefix)mine $(build)$(prefix)optimize $(build)$(prefix)predict $(build)$(prefix)prioritize $(build)$(prefix)similarity $(build)$(prefix)subgraphs $(build)$(prefix)validate

all: fragments $(build)$(prefix)aminoacids $(build)inspire

$(build)$(prefix)classify $(build)$(prefix)exclude $(build)$(prefix)optimize $(build)$(prefix)predict $(build)$(prefix)prioritize $(build)$(prefix)subgraphs: $(build)$(prefix)%: backend/octree.h backend/features.h
$(build)$(prefix)classify $(build)$(prefix)exclude $(build)$(prefix)optimize $(build)$(prefix)predict $(build)$(prefix)prioritize $(build)$(prefix)subgraphs $(build)$(prefix)assign $(build)$(prefix)similarity: $(build)$(prefix)%: frontend/%.cpp backend/%.h common/filesystem.h common/string.h common/exception.h backend/index.h backend/iterators.h backend/protein.h backend/parser.h backend/pdb.h backend/filters.h backend/mmcif.h backend/xml.h common/xml.h
	$(CXX) $(CXXFLAGS) -o $@ $<
	cp ../man/$(patsubst frontend/%.cpp,%,$<).man $(build)$@.1

$(build)$(prefix)fingerprints: $(build)$(prefix)%: backend/subgraphs.h common/graph.h
$(build)$(prefix)fingerprints $(build)$(prefix)filter $(build)$(prefix)validate: $(build)$(prefix)%: backend/index.h backend/features.h backend/octree.h
$(build)$(prefix)fingerprints $(build)$(prefix)filter $(build)$(prefix)validate $(build)$(prefix)index: $(build)$(prefix)%: backend/iterators.h backend/protein.h backend/parser.h backend/pdb.h backend/filters.h backend/mmcif.h backend/xml.h common/xml.h
$(build)$(prefix)fingerprints $(build)$(prefix)filter $(build)$(prefix)validate $(build)$(prefix)index $(build)$(prefix)merge: $(build)$(prefix)%: backend/%.h common/exception.h
$(build)$(prefix)fingerprints $(build)$(prefix)filter $(build)$(prefix)validate $(build)$(prefix)index $(build)$(prefix)merge $(build)$(prefix)aminoacids: $(build)$(prefix)%: frontend/%.cpp common/string.h common/filesystem.h
	$(CXX) $(CXXFLAGS) -o $@ $< -lboost_filesystem
	cp ../man/$(patsubst frontend/%.cpp,%,$<).man $(build)$@.1

$(build)$(prefix)mine: frontend/mine.cpp common/filesystem.h common/string.h common/exception.h backend/mine.h common/multithread.h
	$(CXX) $(CXXFLAGS) -pthread -o $@ $< -lboost_filesystem
	cp ../man/$(patsubst frontend/%.cpp,%,$<).man $(build)$@.1

$(build)$(prefix)features: frontend/features.cpp backend/features.h common/filesystem.h common/string.h backend/index.h backend/iterators.h common/exception.h backend/protein.h backend/parser.h backend/pdb.h backend/filters.h backend/mmcif.h backend/xml.h common/xml.h backend/octree.h common/sasa.h
	$(CXX) $(CXXFLAGS) -pthread -o $@ $< -lboost_filesystem $(rasa)
	cp ../man/$(patsubst frontend/%.cpp,%,$<).man $(build)$@.1

$(build)inspire: frontend/inspire.cpp backend/index.h backend/iterators.h common/exception.h backend/protein.h backend/parser.h common/string.h backend/pdb.h backend/filters.h backend/mmcif.h backend/xml.h common/xml.h common/filesystem.h backend/features.h backend/octree.h backend/subgraphs.h backend/fingerprints.h common/graph.h backend/mine.h common/multithread.h backend/classify.h backend/predict.h common/sasa.h
	$(CXX) $(CXXFLAGS) -pthread -o $@ $< -lboost_filesystem $(rasa)
	cp ../man/$@.man $(build)$@.1

install:
	install -m755 $(wildcard $(build)$(prefix)assign $(build)$(prefix)classify $(build)$(prefix)exclude $(build)$(prefix)features $(build)$(prefix)filter $(build)$(prefix)fingerprints $(build)$(prefix)index $(build)$(prefix)merge $(build)$(prefix)mine $(build)$(prefix)optimize $(build)$(prefix)predict $(build)$(prefix)prioritize $(build)$(prefix)similarity $(build)$(prefix)subgraphs $(build)$(prefix)validate $(build)$(prefix)aminoacids $(build)inspire) $(BINARIES)
	install -m644 $(wildcard $(build)*.1) $(MANPAGES)

uninstall:
	rm -f $(wildcard $(BINARIES)$(prefix)assign $(BINARIES)$(prefix)classify $(BINARIES)$(prefix)exclude $(BINARIES)$(prefix)features $(BINARIES)$(prefix)filter $(BINARIES)$(prefix)fingerprints $(BINARIES)$(prefix)index $(BINARIES)$(prefix)merge $(BINARIES)$(prefix)mine $(BINARIES)$(prefix)optimize $(BINARIES)$(prefix)predict $(BINARIES)$(prefix)prioritize $(BINARIES)$(prefix)similarity $(BINARIES)$(prefix)subgraphs $(BINARIES)$(prefix)validate $(BINARIES)$(prefix)aminoacids $(BINARIES)inspire)
	rm -f $(wildcard $(MANPAGES)$(prefix)assign.1 $(MANPAGES)$(prefix)classify.1 $(MANPAGES)$(prefix)exclude.1 $(MANPAGES)$(prefix)features.1 $(MANPAGES)$(prefix)filter.1 $(MANPAGES)$(prefix)fingerprints.1 $(MANPAGES)$(prefix)index.1 $(MANPAGES)$(prefix)merge.1 $(MANPAGES)$(prefix)mine.1 $(MANPAGES)$(prefix)optimize.1 $(MANPAGES)$(prefix)predict.1 $(MANPAGES)$(prefix)prioritize.1 $(MANPAGES)$(prefix)similarity.1 $(MANPAGES)$(prefix)subgraphs.1 $(MANPAGES)$(prefix)validate.1 $(MANPAGES)$(prefix)aminoacids.1 $(MANPAGES)inspire.1)

clean:
	rm -f $(wildcard $(build)$(prefix)assign $(build)$(prefix)classify $(build)$(prefix)exclude $(build)$(prefix)features $(build)$(prefix)filter $(build)$(prefix)fingerprints $(build)$(prefix)index $(build)$(prefix)merge $(build)$(prefix)mine $(build)$(prefix)optimize $(build)$(prefix)predict $(build)$(prefix)prioritize $(build)$(prefix)similarity $(build)$(prefix)subgraphs $(build)$(prefix)validate $(build)$(prefix)aminoacids $(build)inspire)
	rm -f $(wildcard $(build)*.1)
