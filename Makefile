.DEFAULT_GOAL := all
CXX ?= c++
CXXFLAGS ?= -O3 -g
BCFTOOLS ?= bcftools
PYTHON ?= python3
BUILD_DIR ?= build
PREFIX ?= /usr/local
plugindir ?= $(PREFIX)/libexec/bcftools
EXPORT_DIR ?= build/public-plugin
HEADER_ROOTS = vendor/htslib-1.25 ../locusfinder/vendor/htslib-1.25 /usr/local/include /usr/include
PLUGIN_HTSLIB ?= $(firstword $(foreach root,$(HEADER_ROOTS),$(if $(wildcard $(root)/htslib/hts.h),$(root))))
BCFTOOLS_VERSION ?= $(shell $(BCFTOOLS) --version-only 2>/dev/null | cut -d+ -f1)
PLUGIN_HTSLIB_VERSION ?= $(shell if [ -f "$(PLUGIN_HTSLIB)/version.sh" ]; then cd "$(PLUGIN_HTSLIB)" && sh version.sh; else pkg-config --modversion htslib 2>/dev/null; fi)
PLATFORM ?= $(shell uname -s)
PLUGIN = $(BUILD_DIR)/genecontext.so
HEADERS := $(wildcard src/*.hpp)

ifeq ($(PLATFORM),Darwin)
BCFTOOLS_PATH := $(shell command -v $(BCFTOOLS))
PLUGIN_LDFLAGS ?= -bundle -bundle_loader "$(BCFTOOLS_PATH)" -Wl,-undefined,dynamic_lookup
else
PLUGIN_LDFLAGS ?= -shared
endif

.PHONY: all plugin test install export clean check-deps
all: plugin
plugin: $(PLUGIN)
$(BUILD_DIR):
	mkdir -p "$@"
check-deps:
	@test -f "$(PLUGIN_HTSLIB)/htslib/hts.h" || { echo 'HTSlib headers missing; set PLUGIN_HTSLIB or run bash scripts/bootstrap.sh' >&2; exit 1; }
	@test -n "$(BCFTOOLS_VERSION)" || { echo 'bcftools is required; set BCFTOOLS=/path/to/bcftools' >&2; exit 1; }
	@test -n "$(PLUGIN_HTSLIB_VERSION)" || { echo 'Set PLUGIN_HTSLIB_VERSION to the version of your HTSlib headers' >&2; exit 1; }
$(BUILD_DIR)/plugin.o: src/plugin.cpp $(HEADERS) | $(BUILD_DIR) check-deps
	$(CXX) $(CPPFLAGS) -I"$(PLUGIN_HTSLIB)" $(CXXFLAGS) -std=c++17 -Wall -Wextra -Wpedantic -fPIC -fvisibility=hidden \
	  -DGC_BCFTOOLS_VERSION='"$(BCFTOOLS_VERSION)"' -DGC_HTSLIB_VERSION='"$(PLUGIN_HTSLIB_VERSION)"' -MMD -MP -c $< -o $@
# All HTSlib symbols resolve from the bcftools host: never link a second HTSlib.
$(PLUGIN): $(BUILD_DIR)/plugin.o
	$(CXX) $(CXXFLAGS) $(PLUGIN_LDFLAGS) $(LDFLAGS) $< $(LDLIBS) -o $@
test: plugin
	GC_PLUGIN="$(abspath $(PLUGIN))" BCFTOOLS="$(BCFTOOLS)" $(PYTHON) -B -m unittest discover -s tests -p 'test_*.py' -v
install: plugin
	mkdir -p "$(DESTDIR)$(plugindir)"
	install -m 644 $(PLUGIN) "$(DESTDIR)$(plugindir)/genecontext.so"
export:
	bash scripts/export-plugin.sh "$(EXPORT_DIR)"
clean:
	rm -f $(BUILD_DIR)/plugin.o $(BUILD_DIR)/plugin.d $(PLUGIN)
-include $(BUILD_DIR)/plugin.d
