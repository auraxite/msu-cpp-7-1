CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

run: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

ifeq ($(OS),Windows_NT)
BASH ?= "C:/Program Files/Git/bin/bash.exe"
else
BASH ?= bash
endif

bench:
	$(BASH) bench.sh

ifeq ($(OS),Windows_NT)
clean:
	-del /q run.exe *.bin *.out.tsv 2>nul
else
clean:
	rm -f run *.bin *.out.tsv
endif

.PHONY: bench clean
