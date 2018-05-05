#
# == Paths ==
#
BIN_DIR := bin
SRC_1   := part1
SRC_2   := part2

#
# == Targets ==
#
default: hw2

clean:
	$(RM) $(BUILD_DIR)/*.o $(BIN_DIR)/*

hw2:
	mpicc -o $(BIN_DIR)/$@ -O2 -xHost -qopenmp $(SRC_1)/part1.c
serial:
	mpicc -o $(BIN_DIR)/$@ -O2 -xHost -qopenmp $(SRC_1)/serial.c
