
target = bem_solve

run-parameters = 
debug-mode = on
clean-up-files = *dat

D = $(HOME)/dev/deal.II

# -----------------------------------------------------------------------------
include $D/common/Make.global_options

lib-h-files := $(shell echo $D/include/deal.II/*/*.h)
h-files     := $(wildcard include/*.h)

cc-files    := $(shell echo src/*.cc)
o-files     := $(cc-files:src/%.cc=lib/%.$(OBJEXT))
go-files    := $(cc-files:src/%.cc=lib/%.g.$(OBJEXT))

cc-targets  := $(shell echo *.cc)
o-targets   := $(cc-targets:%.cc=%.$(OBJEXT))
go-targets  := $(cc-targets:%.cc=%.g.$(OBJEXT))

libs.g := $(lib-deal2.g)
libs.o := $(lib-deal2.o)

ifeq ($(debug-mode),on)
  flags     = $(CXXFLAGS.g) -Iinclude
  libraries = $(go-files) $(libs.g)
  targets   = $(go-targets)
  target_o  = lib/$(target).g.$(OBJEXT)
else
  flags     = $(CXXFLAGS.o) -Iinclude
  libraries = $(o-files) $(libs.o)
  targets   = $(o-targets)
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
	$(CXX) -o $@ lib/$(target).g.o $(libraries) $(LIBS) $(LDFLAGS)

run: $(target)$(EXEEXT)
	@echo "============================ Running $<"
	./$(target)$(EXEEXT) $(run-parameters)

clean: clean-lib clean-data
	-rm -f *~ */*~ */*/*~ lib/Makefile.dep

clean-lib:
	-rm -f lib/*.$(OBJEXT) lib/*.g.$(OBJEXT) $(target)$(EXEEXT) lib/TAGS lib/tags

clean-data:
	-rm -f $(clean-up-files)

.PHONY: run clean clean-data clean-lib

lib/Makefile.dep: $(cc-targets) $(cc-files) $(h-files) $(lib-h-files) Makefile
	@echo "============================ Remaking $@"
	$D/common/scripts/make_dependencies $(INCLUDE) -Blib $(cc-targets) $(cc-files) > $@ || (rm -f $@; false)
	if test -s $@; then : else rm -f $@ ; fi

include lib/Makefile.dep
