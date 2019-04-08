BUILD_DIR = $(PWD)/build

default:
	cd $(BUILD_DIR) && make
	cp $(BUILD_DIR)/wave .

clean:
	cd $(BUILD_DIR) && make clean
