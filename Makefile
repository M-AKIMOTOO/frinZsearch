# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS  = -lm -lfftw3f -lstdc++fs #-lfftw3f_threads -lpthread

# Source files
SRCS = \
    src/frinZsearch.cpp   \
    src/frinZargs.cpp     \
    src/frinZread.cpp     \
    src/frinZlogger.cpp   \
    src/frinZparamcal.cpp \
    src/frinZfftshift.cpp \
    src/frinZfitting.cpp

# Object files (derived from SRCS)
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = frinZsearch

# Default target
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)
	@echo ""
	@echo "********************************************************************************"
	@echo "$(TARGET) が正常にビルドされました。"
	@echo ""
	@echo "実行ファイルを PATH/環境変数が通っているディレクトリに移動またはコピーしてください。"
	@echo "例:"
	@echo "  sudo mv/cp $(TARGET) /usr/local/bin/"
	@echo "  mv/cp $(TARGET) ~/bin/  (~/bin がPATHに含まれている場合)"
	@echo ""
	@echo "または、'make install' を実行してインストールすることもできます。"
	@echo "インストール先は Makefile 内の PREFIX や BINDIR 変数で確認・変更できます。"
	@echo "デフォルトのインストール先: $(BINDIR)"
	@echo "********************************************************************************"

# Rule to compile .cpp files to .o files
# This is a pattern rule that applies to all .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# インストールターゲット (オプション)
PREFIX ?= ~
BINDIR ?= $(PREFIX)/bin

install: $(TARGET)
	@echo "Installing $(TARGET) to $(BINDIR)..."
	@mkdir -p $(BINDIR)
	@cp $(TARGET) $(BINDIR)/
	@echo "$(TARGET) が $(BINDIR)/$(TARGET) にインストールされました。"
	@echo "必要に応じて、$(BINDIR) が PATH/環境変数に含まれていることを確認してください。"

uninstall:
	@echo "Uninstalling $(TARGET) from $(BINDIR)..."
	@rm -f $(BINDIR)/$(TARGET)
	@echo "$(TARGET) がアンインストールされました。"

# Target to clean up object files and the executable
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets (targets that are not actual files)
.PHONY: all clean
