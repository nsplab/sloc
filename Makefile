
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
libdeal2-h := $(shell echo $D/include/deal.II/*/*.h)

# -----------------------------------------------------------------------------

# load target source files
include ./Makefile.targets

# stuff from src/
common-o := $(common-cc:src/%.cc=lib/%.$(OBJEXT))
common-go := $(common-cc:src/%.cc=lib/%.g.$(OBJEXT))

# stuff from bin/
bin-o  := $(bin-cc:bin/%.cc=lib/%.$(OBJEXT))
bin-go := $(bin-cc:bin/%.cc=lib/%.g.$(OBJEXT))

# stuff from tests/
tests-o  := $(tests-cc:tests/%.cc=lib/%.$(OBJEXT))
tests-go := $(tests-cc:tests/%.cc=lib/%.g.$(OBJEXT))

# combined source and header files
files-cc := $(common-cc) $(bin-cc) $(tests-cc)
files-h  := $(common-h) $(libdeal2-h)

# complete list of target source files (main target included)
targets-cc := $(bin-cc) $(tests-cc)

# actual targets (strip .cc extension)
targets := $(targets-cc:.cc=$(EXEEXT))


# -----------------------------------------------------------------------------
# Rules

ifeq ($(debug-mode),on)
  flags     = $(CXXFLAGS.g) -Iinclude
  libraries = $(common-go) $(lib-deal2.g)
else
  flags     = $(CXXFLAGS.o) -Iinclude
  libraries = $(common-o) $(lib-deal2.o)
endif

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

lib/Makefile.dep: $(files-cc) $(files-h) Makefile
	@echo "============================ Remaking $@"
	$D/common/scripts/make_dependencies $(INCLUDE) -Blib $(files-cc) > $@ || (rm -f $@; false)
	if test -s $@; then : else rm -f $@ ; fi

include lib/Makefile.dep
