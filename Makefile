
# main target
target = bem_solve

run-parameters = 
debug-mode = on
clean-up-files = *dat

# location of deal.II installation
D = $(HOME)/dev/deal.II

# -----------------------------------------------------------------------------
include $D/common/Make.global_options

# deal.II include files
lib-h-files := $(shell echo $D/include/deal.II/*/*.h)
h-files     := $(wildcard include/*.h)

# deal.II library (debug & optimized)
libs.g := $(lib-deal2.g)
libs.o := $(lib-deal2.o)

# -----------------------------------------------------------------------------

# load target source files
include ./Makefile.targets

# stuff from src/
cc-files    := $(shell echo src/*.cc)
o-files     := $(cc-files:src/%.cc=lib/%.$(OBJEXT))
go-files    := $(cc-files:src/%.cc=lib/%.g.$(OBJEXT))

# stuff from bin/
bin-o  := $(bin-cc:bin/%.cc=lib/%.$(OBJEXT))
bin-go := $(bin-cc:bin/%.cc=lib/%.g.$(OBJEXT))

# stuff from tests/
tests-o  := $(tests-cc:tests/%.cc=lib/%.$(OBJEXT))
tests-go := $(tests-cc:tests/%.cc=lib/%.g.$(OBJEXT))

ifeq ($(debug-mode),on)
  flags     = $(CXXFLAGS.g) -Iinclude
  libraries = $(go-files) $(libs.g)
  target-o  = lib/$(target).g.$(OBJEXT)
else
  flags     = $(CXXFLAGS.o) -Iinclude
  libraries = $(o-files) $(libs.o)
  target-o  = lib/$(target).$(OBJEXT)
endif

# complete list of target source files (main target included)
cc-targets := $(bin-cc) $(tests-cc)

# actual targets (strip .cc extension)
targets := $(cc-targets:.cc=$(EXEEXT))

# -----------------------------------------------------------------------------
# Rules

# by default, build only the main target
default: bin/$(target)$(EXEEXT)

# run the main target
run: bin/$(target)$(EXEEXT)
	@echo "============================ Running $<"
	./bin/$(target)$(EXEEXT) $(run-parameters)

# rule for building everything
all: $(targets)

# -----------------------------------------------------------------------------

# rule for debug object files originating from src/
lib/%.g.$(OBJEXT): src/%.cc
	@echo "==============debug========= $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

# rule for optimized object files originating from src/
lib/%.$(OBJEXT): src/%.cc
	@echo "==============optimized===== $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

# -----------------------------------------------------------------------------

# generic rule for targets in bin/
bin/%: lib/%.g.$(OBJEXT) $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ $< $(libraries) $(LIBS) $(LDFLAGS)

# rule for debug object files originating from bin/
lib/%.g.$(OBJEXT): bin/%.cc
	@echo "==============debug========= $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

# -----------------------------------------------------------------------------

# generic rule for targets in tests/
tests/%: lib/%.g.$(OBJEXT) $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ $< $(libraries) $(LIBS) $(LDFLAGS)

# rule for debug object files originating from tests/
lib/%.g.$(OBJEXT): tests/%.cc
	@echo "==============debug========= $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@


# -----------------------------------------------------------------------------

clean: clean-lib clean-data
	-rm -f *~ */*~ */*/*~ lib/Makefile.dep

clean-lib:
	-rm -f lib/*.$(OBJEXT) lib/*.g.$(OBJEXT) $(target)$(EXEEXT) $(targets) lib/TAGS lib/tags

clean-data:
	-rm -f $(clean-up-files)

.PHONY: run all clean clean-data clean-lib

lib/Makefile.dep: $(cc-targets) $(cc-files) $(h-files) $(lib-h-files) Makefile
	@echo "============================ Remaking $@"
	$D/common/scripts/make_dependencies $(INCLUDE) -Blib $(cc-targets) $(cc-files) > $@ || (rm -f $@; false)
	if test -s $@; then : else rm -f $@ ; fi

include lib/Makefile.dep
