# Alternative GNU Make workspace makefile autogenerated by Premake

ifndef config
  config=debug_x64
endif

ifndef verbose
  SILENT = @
endif

ifeq ($(config),debug_x64)
  lbmApp_config = debug_x64

else ifeq ($(config),release_x64)
  lbmApp_config = release_x64

else ifeq ($(config),relwithdebinfo_x64)
  lbmApp_config = relwithdebinfo_x64

else
  $(error "invalid configuration $(config)")
endif

PROJECTS := lbmApp

.PHONY: all clean help $(PROJECTS) 

all: $(PROJECTS)

lbmApp:
ifneq (,$(lbmApp_config))
	@echo "==== Building lbmApp ($(lbmApp_config)) ===="
	@${MAKE} --no-print-directory -C build -f Makefile config=$(lbmApp_config)
endif

clean:
	@${MAKE} --no-print-directory -C build -f Makefile clean

help:
	@echo "Usage: make [config=name] [target]"
	@echo ""
	@echo "CONFIGURATIONS:"
	@echo "  debug_x64"
	@echo "  release_x64"
	@echo "  relwithdebinfo_x64"
	@echo ""
	@echo "TARGETS:"
	@echo "   all (default)"
	@echo "   clean"
	@echo "   lbmApp"
	@echo ""
	@echo "For more information, see https://github.com/premake/premake-core/wiki"