
# main target
target = bem_solve

run-parameters = 
debug-mode = on
clean-up-files = *dat

# location of deal.II installation
D = $(HOME)/dev/deal.II

# -----------------------------------------------------------------------------
include $D/common/Make.global_options

lib-h-files := $(shell echo $D/include/deal.II/*/*.h)
h-files     := $(wildcard include/*.h)

cc-files    := $(shell echo src/*.cc)
o-files     := $(cc-files:src/%.cc=lib/%.$(OBJEXT))
go-files    := $(cc-files:src/%.cc=lib/%.g.$(OBJEXT))

libs.g := $(lib-deal2.g)
libs.o := $(lib-deal2.o)

ifeq ($(debug-mode),on)
  flags     = $(CXXFLAGS.g) -Iinclude
  libraries = $(go-files) $(libs.g)
  target_o  = lib/$(target).g.$(OBJEXT)
else
  flags     = $(CXXFLAGS.o) -Iinclude
  libraries = $(o-files) $(libs.o)
  target_o  = lib/$(target).$(OBJEXT)
endif

lib/%.g.$(OBJEXT):
	@echo "==============debug========= $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

lib/%.$(OBJEXT):
	@echo "==============optimized===== $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

$(target)$(EXEEXT): $(target_o) $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ $(target_o) $(libraries) $(LIBS) $(LDFLAGS)

run: $(target)$(EXEEXT)
	@echo "============================ Running $<"
	./$(target)$(EXEEXT) $(run-parameters)

# -----------------------------------------------------------------------------
# Other targets

cc-targets := $(shell echo *.cc)
targets    := $(cc-targets:.cc=)

fem_solve: lib/fem_solve.g.o $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ lib/fem_solve.g.o $(libraries) $(LIBS) $(LDFLAGS)

sphere: lib/sphere.g.o $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ lib/sphere.g.o $(libraries) $(LIBS) $(LDFLAGS)

skull: lib/skull.g.o $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ lib/skull.g.o $(libraries) $(LIBS) $(LDFLAGS)

all: $(targets)

# -----------------------------------------------------------------------------

clean: clean-lib clean-data
	-rm -f *~ */*~ */*/*~ lib/Makefile.dep

clean-lib:
	-rm -f lib/*.$(OBJEXT) lib/*.g.$(OBJEXT) $(target)$(EXEEXT) lib/TAGS lib/tags

clean-data:
	-rm -f $(clean-up-files)

.PHONY: run all clean clean-data clean-lib

lib/Makefile.dep: $(cc-targets) $(cc-files) $(h-files) $(lib-h-files) Makefile
	@echo "============================ Remaking $@"
	$D/common/scripts/make_dependencies $(INCLUDE) -Blib $(cc-targets) $(cc-files) > $@ || (rm -f $@; false)
	if test -s $@; then : else rm -f $@ ; fi

include lib/Makefile.dep
