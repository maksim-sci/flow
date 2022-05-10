.PHONY: conf
conf:
	cmake -S . -B build

.PHONY: build
build:
	cmake --build build

.PHONY: test
test:
	cmake --build build --target test

.PHONY: run
run:
	./$(pwd)/build/apps/app