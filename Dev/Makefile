# SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Standalone build for the SoA particle prototype (no Chombo required).

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

TESTS = test_ParticleSoA test_ParticleLayoutMacro

.PHONY: all check clean

all: $(TESTS)

test_ParticleSoA: test_ParticleSoA.cpp CD_ParticleSoA.H CD_RealProxy.H
	$(CXX) $(CXXFLAGS) $< -o $@

test_ParticleLayoutMacro: test_ParticleLayoutMacro.cpp CD_ParticleLayoutMacro.H CD_ParticleSoA.H CD_RealProxy.H
	$(CXX) $(CXXFLAGS) $< -o $@

check: $(TESTS)
	./test_ParticleSoA
	./test_ParticleLayoutMacro

clean:
	$(RM) $(TESTS)
