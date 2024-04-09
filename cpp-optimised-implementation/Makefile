
.PHONY: all
all: main

SOURCES = main.cpp
HEADERS = montgomery.hpp fp2.hpp isog.hpp scallop.hpp

main: $(SOURCES) $(HEADERS)
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -march=native -DNDEBUG \
		$(SOURCES) \
		-lntl \
		-lgmp -lm \
		-lfplll -lmpfr \
		-o $@

debug: $(SOURCES) $(HEADERS)
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -g \
		$(SOURCES) \
		-lntl \
		-lgmp -lm \
		-lfplll -lmpfr \
		-o $@

.PHONY: clean
clean:
	rm -f main debug