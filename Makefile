# ---------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
CC  = g++
CFLAGS    = -g -MMD -MP -Wall -Wextra -Winit-self -Wno-missing-field-initializers
ifeq "$(shell getconf LONG_BIT)" "64"
  LDFLAGS =
else
  LDFLAGS =
endif
LIBS      =
INCLUDE   = -I./include
#TARGET    = ./bin/$(shell basename `readlink -f .`)
TARGETDIR = ./bin
TARGET = $(TARGETDIR)/lbm-sim
SRCDIR    = ./src
ifeq "$(strip $(SRCDIR))" ""
  SRCDIR  = .
endif
SOURCES   = $(wildcard $(SRCDIR)/*.cpp)
OBJDIR    = ./obj
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR  = .
endif
OBJECTS   = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
DEPENDS   = $(OBJECTS:.o=.d)
# -------------------------------------------------------------------------------------------------

all: $(TARGET)

init:
	-mkdir -p img


$(TARGET): $(OBJECTS) $(LIBS)
	-mkdir -p $(TARGETDIR)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<


.PHONY: clean
clean:
	-rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

-include $(DEPENDS)
