
target = bem_solve
run-parameters = 
debug-mode = on
clean-up-files = *dat

D = $(HOME)/dev/deal.II

# -----------------------------------------------------------------------------
include $D/common/Make.global_options

cc-files    := $(shell echo *.cc)
o-files     := $(cc-files:%.cc=%.$(OBJEXT))
go-files    := $(cc-files:%.cc=%.g.$(OBJEXT))
h-files     := $(wildcard *.h)
lib-h-files := $(shell echo $D/include/deal.II/*/*.h)

libs.g := $(lib-deal2.g)
libs.o := $(lib-deal2.o)

ifeq ($(debug-mode),on)
  libraries = $(go-files) $(libs.g)
  flags     = $(CXXFLAGS.g) -I.
else
  libraries = $(o-files) $(libs.o)
  flags     = $(CXXFLAGS.o) -I.
endif

./%.g.$(OBJEXT):
	@echo "==============debug========= $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

./%.$(OBJEXT):
	@echo "==============optimized===== $(<F)  ->  $@"
	$(CXX) $(flags) -c $< -o $@

$(target)$(EXEEXT): $(libraries) Makefile
	@echo "============================ Linking $@"
	$(CXX) -o $@ $(libraries) $(LIBS) $(LDFLAGS)

run: $(target)$(EXEEXT)
	@echo "============================ Running $<"
	./$(target)$(EXEEXT) $(run-parameters)

clean: clean-lib clean-data
	-rm -f *~ */*~ */*/*~ Makefile.dep

clean-lib:
	-rm -f *.$(OBJEXT) *.g.$(OBJEXT) $(target)$(EXEEXT) tags TAGS

clean-data:
	-rm -f $(clean-up-files)

.PHONY: run clean clean-data clean-lib

Makefile.dep: $(target).cc Makefile $(shell echo $D/include/deal.II/*/*.h)
	@echo "============================ Remaking $@"
	$D/common/scripts/make_dependencies $(INCLUDE) -B. $(cc-files) > $@ || (rm -f $@; false)
	if test -s $@; then : else rm -f $@ ; fi

include Makefile.dep
