# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS  = -lm -lfftw3f #-lstdc++fs #-lfftw3f_threads -lpthread

# Executable names
TARGET_SEARCH = frinZsearch
TARGET_RAWVIS = frinZrawvis
TARGETS = $(TARGET_SEARCH) $(TARGET_RAWVIS)

# Source files for frinZsearch
SRCS_SEARCH = \
    src/cpp/frinZsearch.cpp   \
    src/cpp/frinZargs.cpp     \
    src/cpp/frinZread.cpp     \
    src/cpp/frinZlogger.cpp   \
    src/cpp/frinZparamcal.cpp \
    src/cpp/frinZfftshift.cpp \
    src/cpp/frinZfitting.cpp

# Source files for frinZread
SRCS_RAWVIS = \
    src/cpp/frinZrawvis.cpp   \
    src/cpp/frinZread.cpp     \
    src/cpp/frinZoutput.cpp   \
    src/cpp/frinZlogger.cpp   \
    src/cpp/frinZargs.cpp

# All unique object files
OBJS_SEARCH = $(SRCS_SEARCH:.cpp=.o)
OBJS_RAWVIS = $(SRCS_RAWVIS:.cpp=.o)
OBJS = $(sort $(OBJS_SEARCH) $(OBJS_RAWVIS))

# Default target
all: $(TARGETS)

# Rule to link frinZsearch
$(TARGET_SEARCH): $(OBJS_SEARCH)
	$(CXX) $(OBJS_SEARCH) -o $(TARGET_SEARCH) $(LDFLAGS)
	@echo ""
	@echo "********************************************************************************"
	@echo "$(TARGET_SEARCH) が正常にビルドされました。"
	@echo "********************************************************************************"

# Rule to link frinZread
$(TARGET_RAWVIS): $(OBJS_RAWVIS)
	$(CXX) $(OBJS_RAWVIS) -o $(TARGET_RAWVIS) $(LDFLAGS)
	@echo ""
	@echo "********************************************************************************"
	@echo "$(TARGET_RAWVIS) が正常にビルドされました。"
	@echo "********************************************************************************"

# Rule to compile .cpp files to .o files
src/cpp/%.o: src/cpp/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# インストールターゲット (オプション)
PREFIX ?= ~
BINDIR ?= $(PREFIX)/bin

install: $(TARGETS)
	@echo "Installing $(TARGETS) to $(BINDIR)..."
	@mkdir -p $(BINDIR)
	@for T in $(TARGETS); do \
		cp $$T $(BINDIR)/; \
	done
	@echo "$(TARGETS) が $(BINDIR)/ にインストールされました。"
	@echo "必要に応じて、$(BINDIR) が PATH/環境変数に含まれていることを確認してください。"

uninstall:
	@echo "Uninstalling $(TARGETS) from $(BINDIR)..."
	@for T in $(TARGETS); do \
		rm -f $(BINDIR)/$$T; \
	done
	@echo "$(TARGETS) がアンインストールされました。"

# Target to clean up object files and the executable
clean:
	rm -f $(OBJS) $(TARGETS)

# Phony targets (targets that are not actual files)
.PHONY: all clean install uninstall
