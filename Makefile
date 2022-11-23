TARGET_EXEC := pef

TARGET_FILE := ./src/main.cpp
BUILD_DIR := ./build
STD_VERSION := -std=c++11
CXX := g++
SDSL_INCLUDE := ./external/cpp_libs/include
SDSL_LIB := ./external/cpp_libs/lib
LIBS := -I${SDSL_INCLUDE} -I./src/include -L${SDSL_LIB}

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
#CPPFLAGS := -fmodules-ts $(INC_FLAGS) -MMD -MP -Wall -Werror
CPPFLAGS := -MMD -MP -Wall -Weffc++ -Wextra -Wsign-conversion


# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(TARGET_FILE)
	$(CXX) $(STD_VERSION) $(LIBS) $(TARGET_FILE) -o $@ -O3 -lsdsl -ldivsufsort -ldivsufsort64

$(TARGET_FILE):
	mkdir build


.PHONY: clean
clean:
	rm -r $(BUILD_DIR)

cache-clean:
	rm -r gcm.cache

